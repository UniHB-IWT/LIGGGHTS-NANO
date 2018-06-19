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
    This file is from LAMMPS
    LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
    http://lammps.sandia.gov, Sandia National Laboratories
    Steve Plimpton, sjplimp@sandia.gov

    Copyright (2003) Sandia Corporation.  Under the terms of Contract
    DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
    certain rights in this software.  This software is distributed under
    the GNU General Public License.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Naveen Michaud-Agrawal (Johns Hopkins U)
     K-space terms added by Stan Moore (BYU)
------------------------------------------------------------------------- */

#include <mpi.h>
#include <string.h>
#include "compute_group_group.h"
#include "atom.h"
#include "update.h"
#include "force.h"
#include "pair.h"
#include "neighbor.h"
#include "neigh_request.h"
#include "neigh_list.h"
#include "group.h"
#include "kspace.h"
#include "error.h"
#include <cmath>
#include "comm.h"
#include "domain.h"
#include "math_const.h"

using namespace LAMMPS_NS;
using namespace MathConst;

#define SMALL 0.00001

/* ---------------------------------------------------------------------- */

ComputeGroupGroup::ComputeGroupGroup(LAMMPS *lmp, int &iarg, int narg, char **arg) :
  Compute(lmp, iarg, narg, arg)
{
  if (narg < iarg+1) error->all(FLERR,"Illegal compute group/group command");

  scalar_flag = vector_flag = 1;
  size_vector = 3;
  extscalar = 1;
  extvector = 1;

  int n = strlen(arg[iarg]) + 1;
  group2 = new char[n];
  strcpy(group2,arg[iarg]);
  iarg++;

  jgroup = group->find(group2);
  if (jgroup == -1)
    error->all(FLERR,"Compute group/group group ID does not exist");
  jgroupbit = group->bitmask[jgroup];

  pairflag = 1;
  kspaceflag = 0;
  boundaryflag = 1;

  while (iarg < narg) {
    if (strcmp(arg[iarg],"pair") == 0) {
      if (iarg+2 > narg)
        error->all(FLERR,"Illegal compute group/group command");
      if (strcmp(arg[iarg+1],"yes") == 0) pairflag = 1;
      else if (strcmp(arg[iarg+1],"no") == 0) pairflag = 0;
      else error->all(FLERR,"Illegal compute group/group command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"kspace") == 0) {
      if (iarg+2 > narg)
        error->all(FLERR,"Illegal compute group/group command");
      if (strcmp(arg[iarg+1],"yes") == 0) kspaceflag = 1;
      else if (strcmp(arg[iarg+1],"no") == 0) kspaceflag = 0;
      else error->all(FLERR,"Illegal compute group/group command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"boundary") == 0) {
      if (iarg+2 > narg)
        error->all(FLERR,"Illegal compute group/group command");
      if (strcmp(arg[iarg+1],"yes") == 0) boundaryflag = 1;
      else if (strcmp(arg[iarg+1],"no") == 0) boundaryflag  = 0;
      else error->all(FLERR,"Illegal compute group/group command");
      iarg += 2;
    } else error->all(FLERR,"Illegal compute group/group command");
  }

  vector = new double[3];
}

/* ---------------------------------------------------------------------- */

ComputeGroupGroup::~ComputeGroupGroup()
{
  delete [] group2;
  delete [] vector;
}

/* ---------------------------------------------------------------------- */

void ComputeGroupGroup::init()
{
  // if non-hybrid, then error if single_enable = 0
  // if hybrid, let hybrid determine if sub-style sets single_enable = 0

  if (pairflag && force->pair == NULL)
    error->all(FLERR,"No pair style defined for compute group/group");
  if (force->pair_match("hybrid",0) == NULL && force->pair->single_enable == 0)
    error->all(FLERR,"Pair style does not support compute group/group");

  // error if Kspace style does not compute group/group interactions

  if (kspaceflag && force->kspace == NULL)
    error->all(FLERR,"No Kspace style defined for compute group/group");
  if (kspaceflag && force->kspace->group_group_enable == 0)
    error->all(FLERR,"Kspace style does not support compute group/group");

  if (pairflag) {
    pair = force->pair;
    cutsq = force->pair->cutsq;
  } else pair = NULL;

  if (kspaceflag) kspace = force->kspace;
  else kspace = NULL;

  // compute Kspace correction terms

  if (kspaceflag) {
    kspace_correction();
    if (fabs(e_correction) > SMALL && comm->me == 0) {
      char str[512];
      sprintf(str,"Both groups in compute group/group have a net charge; "
              "the Kspace boundary correction to energy will be non-zero");
      error->warning(FLERR,str);
    }
  }

  // recheck that group 2 has not been deleted

  jgroup = group->find(group2);
  if (jgroup == -1)
    error->all(FLERR,"Compute group/group group ID does not exist");
  jgroupbit = group->bitmask[jgroup];

  // need an occasional half neighbor list

  if (pairflag) {
    int irequest = neighbor->request((void *) this);
    neighbor->requests[irequest]->pair = 0;
    neighbor->requests[irequest]->compute = 1;
    neighbor->requests[irequest]->occasional = 1;
  }
}

/* ---------------------------------------------------------------------- */

void ComputeGroupGroup::init_list(int id, NeighList *ptr)
{
  list = ptr;
}

/* ---------------------------------------------------------------------- */

double ComputeGroupGroup::compute_scalar()
{
  invoked_scalar = invoked_vector = update->ntimestep;

  scalar = 0.0;
  vector[0] = vector[1] = vector[2] = 0.0;

  if (pairflag) pair_contribution();
  if (kspaceflag) kspace_contribution();

  return scalar;
}

/* ---------------------------------------------------------------------- */

void ComputeGroupGroup::compute_vector()
{
  invoked_scalar = invoked_vector = update->ntimestep;

  scalar = 0.0;
  vector[0] = vector[1] = vector[2] = 0.0;

  if (pairflag) pair_contribution();
  if (kspaceflag) kspace_contribution();
}

/* ---------------------------------------------------------------------- */

void ComputeGroupGroup::pair_contribution()
{
  int i,j,ii,jj,inum,jnum,itype,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz;
  double rsq,eng,fpair,factor_coul,factor_lj;
  int *ilist,*jlist,*numneigh,**firstneigh;

  double **x = atom->x;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  double *special_coul = force->special_coul;
  double *special_lj = force->special_lj;
  int newton_pair = force->newton_pair;

  // invoke half neighbor list (will copy or build if necessary)

  neighbor->build_one(list->index);

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // loop over neighbors of my atoms
  // skip if I,J are not in 2 groups

  double one[4];
  one[0] = one[1] = one[2] = one[3] = 0.0;

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    if (!(mask[i] & groupbit || mask[i] & jgroupbit)) continue; // skip if atom I is not in either group

    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      factor_lj = special_lj[sbmask(j)];
      factor_coul = special_coul[sbmask(j)];
      j &= NEIGHMASK;

      if (!(mask[j] & groupbit || mask[j] & jgroupbit)) continue; // skip if atom J is not in either group

      int ij_flag = 0;
      int ji_flag = 0;
      if (mask[i] & groupbit && mask[j] & jgroupbit) ij_flag = 1;
      if (mask[j] & groupbit && mask[i] & jgroupbit) ji_flag = 1;
      if (!ij_flag && !ji_flag) continue; // skip if atoms I,J are only in the same group

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      jtype = type[j];

      if (rsq < cutsq[itype][jtype]) {
        eng = pair->single(i,j,itype,jtype,rsq,factor_coul,factor_lj,fpair);

        // energy only computed once so tally full amount
        // force tally is jgroup acting on igroup

        if (newton_pair || j < nlocal) {
          one[0] += eng;
          if (ij_flag) {
            one[1] += delx*fpair;
            one[2] += dely*fpair;
            one[3] += delz*fpair;
          }
          if (ji_flag) {
            one[1] -= delx*fpair;
            one[2] -= dely*fpair;
            one[3] -= delz*fpair;
          }

        // energy computed twice so tally half amount
        // only tally force if I own igroup atom

        } else {
          one[0] += 0.5*eng;
          if (ij_flag) {
            one[1] += delx*fpair;
            one[2] += dely*fpair;
            one[3] += delz*fpair;
          }
        }
      }
    }
  }

  double all[4];
  MPI_Allreduce(one,all,4,MPI_DOUBLE,MPI_SUM,world);
  scalar += all[0];
  vector[0] += all[1]; vector[1] += all[2]; vector[2] += all[3];
}

/* ---------------------------------------------------------------------- */

void ComputeGroupGroup::kspace_contribution()
{
  double *vector_kspace = force->kspace->f2group;

  force->kspace->compute_group_group(groupbit,jgroupbit,0);
  scalar += 2.0*force->kspace->e2group;
  vector[0] += vector_kspace[0];
  vector[1] += vector_kspace[1];
  vector[2] += vector_kspace[2];

  // subtract extra A <--> A Kspace interaction so energy matches
  //   real-space style of compute group-group
  // add extra Kspace term to energy

  force->kspace->compute_group_group(groupbit,jgroupbit,1);
  scalar -= force->kspace->e2group;

  // self energy correction term

  scalar -= e_self;

  // k=0 boundary correction term

  if (boundaryflag) {
    double xprd = domain->xprd;
    double yprd = domain->yprd;
    double zprd = domain->zprd;

    // adjustment of z dimension for 2d slab Ewald
    // 3d Ewald just uses zprd since slab_volfactor = 1.0

    double volume = xprd*yprd*zprd*force->kspace->slab_volfactor;
    scalar -= e_correction/volume;
  }
}

/* ---------------------------------------------------------------------- */

void ComputeGroupGroup::kspace_correction()
{

  // total charge of groups A & B, needed for correction term

  double qsqsum_group,qsum_A,qsum_B;
  qsqsum_group = qsum_A = qsum_B = 0.0;

  double *q = atom->q;
  int *mask = atom->mask;
  int groupbit_A = groupbit;
  int groupbit_B = jgroupbit;

  for (int i = 0; i < atom->nlocal; i++) {
    if ((mask[i] & groupbit_A) && (mask[i] & groupbit_B))
      qsqsum_group += q[i]*q[i];
    if (mask[i] & groupbit_A) qsum_A += q[i];
    if (mask[i] & groupbit_B) qsum_B += q[i];
  }

  double tmp;
  MPI_Allreduce(&qsqsum_group,&tmp,1,MPI_DOUBLE,MPI_SUM,world);
  qsqsum_group = tmp;

  MPI_Allreduce(&qsum_A,&tmp,1,MPI_DOUBLE,MPI_SUM,world);
  qsum_A = tmp;

  MPI_Allreduce(&qsum_B,&tmp,1,MPI_DOUBLE,MPI_SUM,world);
  qsum_B = tmp;

  double g_ewald = force->kspace->g_ewald;

  double scale = 1.0;
  const double qscale = force->qqrd2e * scale;

  // self-energy correction

  e_self = qscale * g_ewald*qsqsum_group/MY_PIS;
  e_correction = 2.0*qsum_A*qsum_B;

  // subtract extra AA terms

  qsum_A = qsum_B = 0.0;

  for (int i = 0; i < atom->nlocal; i++) {
    if (!((mask[i] & groupbit_A) && (mask[i] & groupbit_B)))
      continue;

    if (mask[i] & groupbit_A) qsum_A += q[i];
    if (mask[i] & groupbit_B) qsum_B += q[i];
  }

  MPI_Allreduce(&qsum_A,&tmp,1,MPI_DOUBLE,MPI_SUM,world);
  qsum_A = tmp;

  MPI_Allreduce(&qsum_B,&tmp,1,MPI_DOUBLE,MPI_SUM,world);
  qsum_B = tmp;

  // k=0 energy correction term (still need to divide by volume above)

  e_correction -= qsum_A*qsum_B;
  e_correction *= qscale * MY_PI2 / (g_ewald*g_ewald);
}
