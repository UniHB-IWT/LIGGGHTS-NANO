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

#include <string.h>
#include "compute_msd.h"
#include "atom.h"
#include "update.h"
#include "group.h"
#include "domain.h"
#include "modify.h"
#include "fix_store.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeMSD::ComputeMSD(LAMMPS *lmp, int &iarg, int narg, char **arg) :
  Compute(lmp, iarg, narg, arg)
{
  if (narg < iarg) error->all(FLERR,"Illegal compute msd command");

  vector_flag = 1;
  size_vector = 4;
  extvector = 0;

  // optional args

  comflag = 0;

  while (iarg < narg) {
    if (strcmp(arg[iarg],"com") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal compute msd command");
      if (strcmp(arg[iarg+1],"no") == 0) comflag = 0;
      else if (strcmp(arg[iarg+1],"yes") == 0) comflag = 1;
      else error->all(FLERR,"Illegal compute msd command");
      iarg += 2;
    } else error->all(FLERR,"Illegal compute msd command");
  }

  // create a new fix STORE style
  // id = compute-ID + COMPUTE_STORE, fix group = compute group

  int n = strlen(id) + strlen("_COMPUTE_STORE") + 1;
  id_fix = new char[n];
  strcpy(id_fix,id);
  strcat(id_fix,"_COMPUTE_STORE");

  char **newarg = new char*[5];
  newarg[0] = id_fix;
  newarg[1] = group->names[igroup];
  newarg[2] = (char *) "STORE";
  newarg[3] = (char *) "1";
  newarg[4] = (char *) "3";
  modify->add_fix(5,newarg);
  fix = (FixStore *) modify->fix[modify->nfix-1];
  delete [] newarg;

  // calculate xu,yu,zu for fix store array
  // skip if reset from restart file

  if (fix->restart_reset) fix->restart_reset = 0;
  else {
    double **xoriginal = fix->astore;

    double **x = atom->x;
    int *mask = atom->mask;
    tagint *image = atom->image;
    int nlocal = atom->nlocal;

    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) domain->unmap(x[i],image[i],xoriginal[i]);
      else xoriginal[i][0] = xoriginal[i][1] = xoriginal[i][2] = 0.0;

    // adjust for COM if requested

    if (comflag) {
      double cm[3];
      masstotal = group->mass(igroup);
      group->xcm(igroup,masstotal,cm);
      for (int i = 0; i < nlocal; i++)
        if (mask[i] & groupbit) {
          xoriginal[i][0] -= cm[0];
          xoriginal[i][1] -= cm[1];
          xoriginal[i][2] -= cm[2];
        }
    }
  }

  // displacement vector

  vector = new double[4];
}

/* ---------------------------------------------------------------------- */

ComputeMSD::~ComputeMSD()
{
  // check nfix in case all fixes have already been deleted

  if (modify->nfix) modify->delete_fix(id_fix);

  delete [] id_fix;
  delete [] vector;
}

/* ---------------------------------------------------------------------- */

void ComputeMSD::init()
{
  // set fix which stores original atom coords

  int ifix = modify->find_fix(id_fix);
  if (ifix < 0) error->all(FLERR,"Could not find compute msd fix ID");
  fix = (FixStore *) modify->fix[ifix];

  // nmsd = # of atoms in group

  nmsd = group->count(igroup);
  masstotal = group->mass(igroup);
}

/* ---------------------------------------------------------------------- */

void ComputeMSD::compute_vector()
{
  invoked_vector = update->ntimestep;

  // cm = current center of mass

  double cm[3];
  if (comflag) group->xcm(igroup,masstotal,cm);
  else cm[0] = cm[1] = cm[2] = 0.0;

  // dx,dy,dz = displacement of atom from original position
  // original unwrapped position is stored by fix
  // relative to center of mass if comflag is set
  // for triclinic, need to unwrap current atom coord via h matrix

  double **xoriginal = fix->astore;

  double **x = atom->x;
  int *mask = atom->mask;
  tagint *image = atom->image;
  int nlocal = atom->nlocal;

  double *h = domain->h;
  double xprd = domain->xprd;
  double yprd = domain->yprd;
  double zprd = domain->zprd;

  double dx,dy,dz;
  int xbox,ybox,zbox;

  double msd[4];
  msd[0] = msd[1] = msd[2] = msd[3] = 0.0;

  if (domain->triclinic == 0) {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        xbox = (image[i] & IMGMASK) - IMGMAX;
        ybox = (image[i] >> IMGBITS & IMGMASK) - IMGMAX;
        zbox = (image[i] >> IMG2BITS) - IMGMAX;
        dx = x[i][0] + xbox*xprd - cm[0] - xoriginal[i][0];
        dy = x[i][1] + ybox*yprd - cm[1] - xoriginal[i][1];
        dz = x[i][2] + zbox*zprd - cm[2] - xoriginal[i][2];
        msd[0] += dx*dx;
        msd[1] += dy*dy;
        msd[2] += dz*dz;
        msd[3] += dx*dx + dy*dy + dz*dz;
      }

  } else {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        xbox = (image[i] & IMGMASK) - IMGMAX;
        ybox = (image[i] >> IMGBITS & IMGMASK) - IMGMAX;
        zbox = (image[i] >> IMG2BITS) - IMGMAX;
        dx = x[i][0] + h[0]*xbox + h[5]*ybox + h[4]*zbox -
          cm[0] - xoriginal[i][0];
        dy = x[i][1] + h[1]*ybox + h[3]*zbox - cm[1] - xoriginal[i][1];
        dz = x[i][2] + h[2]*zbox - cm[2] - xoriginal[i][2];
        msd[0] += dx*dx;
        msd[1] += dy*dy;
        msd[2] += dz*dz;
        msd[3] += dx*dx + dy*dy + dz*dz;
      }
  }

  MPI_Allreduce(msd,vector,4,MPI_DOUBLE,MPI_SUM,world);
  if (nmsd) {
    vector[0] /= nmsd;
    vector[1] /= nmsd;
    vector[2] /= nmsd;
    vector[3] /= nmsd;
  }
}
