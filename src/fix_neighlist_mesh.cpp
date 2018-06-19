/* ----------------------------------------------------------------------
    This is the

    ██╗     ██╗ ██████╗  ██████╗  ██████╗ ██╗  ██╗████████╗███████╗
    ██║     ██║██╔════╝ ██╔════╝ ██╔════╝ ██║  ██║╚══██╔══╝██╔════╝
    ██║     ██║██║  ███╗██║  ███╗██║  ███╗███████║   ██║   ███████╗
    ██║     ██║██║   ██║██║   ██║██║   ██║██╔══██║   ██║   ╚════██║
    ███████╗██║╚██████╔╝╚██████╔╝╚██████╔╝██║  ██║   ██║   ███████║
    ╚══════╝╚═╝ ╚═════╝  ╚═════╝  ╚═════╝ ╚═╝  ╚═╝   ╚═╝   ╚══════╝®

    DEM simulation engine, released by
    DCS Computing Gmbh, Linz, Austria
    http://www.dcs-computing.com, office@dcs-computing.com

    LIGGGHTS® is part of CFDEM®project:
    http://www.liggghts.com | http://www.cfdem.com

    Core developer and main author:
    Christoph Kloss, christoph.kloss@dcs-computing.com

    LIGGGHTS® is open-source, distributed under the terms of the GNU Public
    License, version 2 or later. It is distributed in the hope that it will
    be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
    of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. You should have
    received a copy of the GNU General Public License along with LIGGGHTS®.
    If not, see http://www.gnu.org/licenses . See also top-level README
    and LICENSE files.

    LIGGGHTS® and CFDEM® are registered trade marks of DCS Computing GmbH,
    the producer of the LIGGGHTS® software and the CFDEM®coupling software
    See http://www.cfdem.com/terms-trademark-policy for details.

-------------------------------------------------------------------------
    Contributing author and copyright for this file:

    Christoph Kloss (DCS Computing GmbH, Linz)
    Christoph Kloss (JKU Linz)
    Richard Berger (JKU Linz)
    Philippe Seil (JKU Linz)

    Copyright 2012-     DCS Computing GmbH, Linz
    Copyright 2009-2012 JKU Linz
------------------------------------------------------------------------- */

#include "fix_neighlist_mesh.h"
#include "fix_mesh_surface.h"
#include "fix_property_atom.h"
#include "modify.h"
#include "container.h"
#include "bounding_box.h"
#include "neighbor.h"
#include "atom.h"
#include "domain.h"
#include "vector_liggghts.h"
#include "update.h"
#include <stdio.h>
#include <algorithm>
#include "atom_vec_ellipsoid.h"

using namespace LAMMPS_NS;
using namespace FixConst;

#define SMALL_DELTA skin/(70.*M_PI)

FixNeighlistMesh::FixNeighlistMesh(LAMMPS *lmp, int narg, char **arg)
: Fix(lmp,narg,arg),
  fix_nneighs_(0),
  fix_nneighs_name_(0),
  buildNeighList(false),
  numAllContacts_(0),
  globalNumAllContacts_(false),
  mbinx(0),
  mbiny(0),
  mbinz(0),
  maxhead(0),
  bins(NULL),
  binhead(NULL),
  skin(0.0),
  distmax(0.0),
  x(NULL),
  r(NULL),
  changingMesh(false),
  changingDomain(false),
  last_bin_update(-1),
  avec(0),
  otherList_(false)
{
    if(!modify->find_fix_id(arg[3]) || !dynamic_cast<FixMeshSurface*>(modify->find_fix_id(arg[3])))
        error->fix_error(FLERR,this,"illegal caller");

    caller_ = static_cast<FixMeshSurface*>(modify->find_fix_id(arg[3]));
    mesh_ = caller_->triMesh();

    if(5 == narg)
    {
        if(0 == strcmp(arg[4],"other_yes"))
            otherList_ = true;
        else if(0 == strcmp(arg[4],"other_no"))
            otherList_ = false;
        else error->fix_error(FLERR,this,"illegal");
        
    }

    groupbit_wall_mesh = groupbit;
}

/* ---------------------------------------------------------------------- */

FixNeighlistMesh::~FixNeighlistMesh()
{
    delete [] fix_nneighs_name_;
    last_bin_update = -1;
}

/* ---------------------------------------------------------------------- */

void FixNeighlistMesh::post_create()
{
    // register
    if(!fix_nneighs_)
    {
        const char* fixarg[9];
        delete [] fix_nneighs_name_;
        fix_nneighs_name_ = new char[strlen(mesh_->mesh_id())+strlen(id)+1+20];

        if(otherList_)
            sprintf(fix_nneighs_name_,"n_neighs_mesh_%s_fix_%s",mesh_->mesh_id(),id);
        
        else
            sprintf(fix_nneighs_name_,"n_neighs_mesh_%s",mesh_->mesh_id());

        fixarg[0]=fix_nneighs_name_;
        fixarg[1]="all";
        fixarg[2]="property/atom";
        fixarg[3]=fix_nneighs_name_;
        fixarg[4]="scalar"; // 1 vector per particle to be registered
        fixarg[5]="yes";    // restart - REQUIRED!
        fixarg[6]="yes";     // communicate ghost
        fixarg[7]="no";     // communicate rev
        fixarg[8]="0.";
        fix_nneighs_ = modify->add_fix_property_atom(9,const_cast<char**>(fixarg),style);

        fix_nneighs_->just_created = false;
        fix_nneighs_->set_internal();
    }
    //check for aspherical
    avec = (AtomVecEllipsoid *) atom->style_match("ellipsoid");
}

/* ---------------------------------------------------------------------- */

void FixNeighlistMesh::initializeNeighlist()
{
    changingMesh = mesh_->isMoving() || mesh_->isDeforming();
    changingDomain = (domain->nonperiodic == 2) || domain->box_change;

    // remove old lists, init new ones
    
    const size_t nall = mesh_->sizeLocal()+mesh_->sizeGhost();

    while(triangles.size() > nall) {
        triangles.pop_back();
    }

    while(triangles.size() < nall) {
        triangles.push_back(TriangleNeighlist());
    }

    for(size_t iTri = 0; iTri < nall; iTri++) {
        TriangleNeighlist & triangle = triangles[iTri];
        triangle.contacts.reserve(std::max(triangle.contacts.capacity(), static_cast<size_t>(128)));
    }
}

/* ---------------------------------------------------------------------- */

void FixNeighlistMesh::setup_pre_force(int foo)
{
    
    pre_neighbor();
    pre_force(0);
}

/* ---------------------------------------------------------------------- */

void FixNeighlistMesh::min_setup_pre_force(int foo)
{
    
    pre_neighbor();
    pre_force(0);
}

/* ---------------------------------------------------------------------- */

void FixNeighlistMesh::pre_delete(bool unfixflag)
{
    if(unfixflag)
        modify->delete_fix(fix_nneighs_->id);
}

/* ---------------------------------------------------------------------- */

int FixNeighlistMesh::setmask()
{
    int mask = 0;
    mask |= MIN_PRE_NEIGHBOR;
    mask |= PRE_NEIGHBOR;
    mask |= MIN_PRE_FORCE;
    mask |= PRE_FORCE;
    mask |= POST_RUN;
    return mask;
}

/* ---------------------------------------------------------------------- */

void FixNeighlistMesh::pre_neighbor()
{
    buildNeighList = true;
}

/* ---------------------------------------------------------------------- */

void FixNeighlistMesh::min_pre_force(int vflag)
{
    pre_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixNeighlistMesh::pre_force(int)
{
    if(!buildNeighList) return;

    changingMesh = mesh_->isMoving() || mesh_->isDeforming();
    changingDomain = (domain->nonperiodic == 2) || domain->box_change;

    buildNeighList = false;
    numAllContacts_ = 0;

    // set num_neigh = 0
    memset(fix_nneighs_->vector_atom, 0, atom->nmax*sizeof(double));

    x = atom->x;
    r = atom->radius;

    if(neighbor->style != 1)
        error->all(FLERR,"Please use style 'bin' in the 'neighbor' command together with triangular walls");

    double rmax = 0.5*(neighbor->cutneighmax - neighbor->skin);
    double prev_skin = skin;
    double prev_distmax = distmax;

    if(changingMesh)
    {
      skin = neighbor->skin;
      
      distmax = neighbor->cutneighmax + SMALL_DELTA;
    }
    else
    {
      skin = 0.5*neighbor->skin;
      
      distmax = neighbor->cutneighmax - rmax + SMALL_DELTA;
    }

    mbinx = neighbor->mbinx;
    mbiny = neighbor->mbiny;
    mbinz = neighbor->mbinz;
    bins = neighbor->bins;
    binhead = neighbor->binhead;
    maxhead = neighbor->maxhead;

    const size_t nall = mesh_->sizeLocal() + mesh_->sizeGhost();

    // update cache if necessary
    if (triangles.size() != nall) {
      initializeNeighlist();
    }

    // update precomputed bins if necessary
    if((skin != prev_skin) || (distmax != prev_distmax) || (neighbor->last_setup_bins_timestep > last_bin_update)) {
      generate_bin_list(nall);
    }

    // manually trigger binning if no pairwise neigh lists exist
    if(0 == neighbor->n_blist() && bins)
        neighbor->bin_atoms();
    else if(!bins)
        error->one(FLERR,"wrong neighbor setting for fix neighlist/mesh");

    for(size_t iTri = 0; iTri < nall; iTri++) {
      TriangleNeighlist & triangle = triangles[iTri];
      handleTriangle(iTri);
      numAllContacts_ += triangle.contacts.size();
    }

    if(globalNumAllContacts_)
        MPI_Sum_Scalar(numAllContacts_,world);

    fix_nneighs_->do_forward_comm();
}

/* ---------------------------------------------------------------------- */

void FixNeighlistMesh::checkBin(AtomVecEllipsoid::Bonus *bonus, std::vector<int>& neighbors, int& nchecked, double contactDistanceFactor, int *mask, int nlocal, int iBin, int iTri, bool haveNonSpherical, int *ellipsoid, double *shape)
{
    int iAtom = binhead[iBin];

    // only handle local atoms and periodic ghosts
    while(iAtom != -1)
    {
      if((iAtom > nlocal) && (!domain->is_periodic_ghost(iAtom)))
      {
          if(bins) iAtom = bins[iAtom];
          else iAtom = -1;

          continue;
      }

      if(! (mask[iAtom] & groupbit_wall_mesh))
      {
          if(bins) iAtom = bins[iAtom];
          else iAtom = -1;
          continue;
      }
      nchecked++;

      if(0) {}
      #ifdef TRI_LINE_ACTIVE_FLAG
      else if(haveNonSpherical) //if non-spherical, check line interaction as well
      {
          double length;
          double cylRadius;
          shape     = bonus[ellipsoid[iAtom]].shape;
          length    = 2.*MathExtraLiggghts::max(shape[0],shape[1],shape[2]);
          cylRadius =    MathExtraLiggghts::min(shape[0],shape[1],shape[2]);
          if( mesh_->resolveTriSegmentNeighbuild(iTri, x[iAtom], length*contactDistanceFactor, cylRadius, skin ) )
          {
            neighbors.push_back(iAtom);
            fix_nneighs_->set_vector_atom_int(iAtom, fix_nneighs_->get_vector_atom_int(iAtom)+1); // num_neigh++
          }
      }
      #endif
      else if(mesh_->resolveTriSphereNeighbuild(iTri,r ? r[iAtom]*contactDistanceFactor : 0. ,x[iAtom],r ? skin : (distmax+skin) ))
      {
        // include iAtom in neighbor list
        neighbors.push_back(iAtom);
        fix_nneighs_->set_vector_atom_int(iAtom, fix_nneighs_->get_vector_atom_int(iAtom)+1); // num_neigh++
        
      }
      if(bins) iAtom = bins[iAtom];
      else iAtom = -1;
    }
}

/* ---------------------------------------------------------------------- */

void FixNeighlistMesh::handleTriangle(int iTri)
{
    TriangleNeighlist & triangle = triangles[iTri];
    std::vector<int> & neighbors = triangle.contacts;
    int & nchecked = triangle.nchecked;
    int *mask = atom->mask;
    int ixMin(0),ixMax(0),iyMin(0),iyMax(0),izMin(0),izMax(0);
    int nlocal = atom->nlocal;
    double contactDistanceFactor = neighbor->contactDistanceFactor;

    int                     *ellipsoid  = atom->ellipsoid;
    AtomVecEllipsoid::Bonus *bonus = 0;
    double *shape = 0;
    bool    haveNonSpherical = false;
    if(ellipsoid)
    {
        haveNonSpherical = true;
        bonus = avec->bonus;
    }

    neighbors.clear();
    nchecked = 0;

    // only do this if I own particles
    if(nlocal)
    {
      if(changingMesh || changingDomain)
      {
        getBinBoundariesForTriangle(iTri,ixMin,ixMax,iyMin,iyMax,izMin,izMax);
    
        for(int ix=ixMin;ix<=ixMax;ix++) {
          for(int iy=iyMin;iy<=iyMax;iy++) {
            for(int iz=izMin;iz<=izMax;iz++) {
              const int iBin = iz*mbiny*mbinx + iy*mbinx + ix;
              if(iBin < 0 || iBin >= maxhead) continue;
              checkBin(bonus, neighbors, nchecked, contactDistanceFactor, mask, nlocal, iBin, iTri, haveNonSpherical, ellipsoid, shape);
            }
          }
        }
      } else {
        const std::vector<int> & triangleBins = triangle.bins;
        const int bincount = triangleBins.size();
        for(int i = 0; i < bincount; i++) {
          const int iBin = triangleBins[i];
          checkBin(bonus, neighbors, nchecked, contactDistanceFactor, mask, nlocal, iBin, iTri, haveNonSpherical, ellipsoid, shape);
        }
      }
    }

}

/* ---------------------------------------------------------------------- */

void FixNeighlistMesh::getBinBoundariesFromBoundingBox(BoundingBox &b,
      int &ixMin,int &ixMax,int &iyMin,int &iyMax,int &izMin,int &izMax)
{
    double delta = distmax;
    double tri_xmin[3] = {b.xLo-delta,b.yLo-delta,b.zLo-delta};
    double tri_xmax[3] = {b.xHi+delta,b.yHi+delta,b.zHi+delta};

    //int binmin = neighbor->coord2bin(tri_xmin,ixMin,iyMin,izMin);
    neighbor->coord2bin(tri_xmin,ixMin,iyMin,izMin);
    
    //int binmax= neighbor->coord2bin(tri_xmax,ixMax,iyMax,izMax);
    neighbor->coord2bin(tri_xmax,ixMax,iyMax,izMax);
    
}

/* ---------------------------------------------------------------------- */

void FixNeighlistMesh::getBinBoundariesForTriangle(int iTri, int &ixMin,int &ixMax,int &iyMin,int &iyMax,int &izMin,int &izMax)
{
  // disable optimization for movingMesh or shrink-wrapped domain
  if(changingMesh || changingDomain) {
    BoundingBox b = mesh_->getElementBoundingBoxOnSubdomain(iTri);
    // extend bbox by cutneighmax and get bin boundaries
    getBinBoundariesFromBoundingBox(b,ixMin,ixMax,iyMin,iyMax,izMin,izMax);
  } else {
    // use cached boundary information
    const BinBoundary & b = triangles[iTri].boundary;
    ixMin = b.xlo;
    ixMax = b.xhi;
    iyMin = b.ylo;
    iyMax = b.yhi;
    izMin = b.zlo;
    izMax = b.zhi;
  }
}

/* ---------------------------------------------------------------------- */

void FixNeighlistMesh::post_run()
{
  last_bin_update = -1; // reset binning for possible next run
}

/* ---------------------------------------------------------------------- */

void FixNeighlistMesh::generate_bin_list(size_t nall)
{
  // precompute triangle bin boundaries
  // disable optimization for changing mesh or domain
  if (!(changingMesh || changingDomain)) {
    double dx = neighbor->binsizex / 2.0;
    double dy = neighbor->binsizey / 2.0;
    double dz = neighbor->binsizez / 2.0;
    double maxdiag = sqrt(dx * dx + dy * dy + dz * dz);

    for (size_t iTri = 0; iTri < nall; iTri++) {
      TriangleNeighlist & triangle = triangles[iTri];
      std::vector<int> & binlist = triangle.bins;
      binlist.clear();

      BinBoundary& bb = triangle.boundary;
      BoundingBox b = mesh_->getElementBoundingBoxOnSubdomain(iTri);

      // extend bbox by cutneighmax and get bin boundaries
      getBinBoundariesFromBoundingBox(b, bb.xlo, bb.xhi, bb.ylo, bb.yhi, bb.zlo, bb.zhi);

      // look at bins and exclude unnecessary ones
      double center[3];
      int total = 0;
      for (int ix = bb.xlo; ix <= bb.xhi; ix++) {
        for (int iy = bb.ylo; iy <= bb.yhi; iy++) {
          for (int iz = bb.zlo; iz <= bb.zhi; iz++) {
            int iBin = iz * mbiny * mbinx + iy * mbinx + ix;
            if (iBin < 0 || iBin >= maxhead)
              continue;

            // determine center of bin (ix, iy, iz)
            neighbor->bin_center(ix, iy, iz, center);

            if (mesh_->resolveTriSphereNeighbuild(iTri, maxdiag, center, distmax + skin))
            {
              binlist.push_back(iBin);
            }
            total++;
          }
        }
      }
      
    }
  }

  last_bin_update = update->ntimestep;
}

int FixNeighlistMesh::getSizeNumContacts()
{
  return mesh_->sizeLocal() + mesh_->sizeGhost();
}
