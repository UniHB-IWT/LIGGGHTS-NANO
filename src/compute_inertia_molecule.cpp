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

#include "compute_inertia_molecule.h"
#include "atom.h"
#include "update.h"
#include "domain.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeInertiaMolecule::ComputeInertiaMolecule(LAMMPS *lmp, int &iarg, int narg, char **arg) :
  Compute(lmp, iarg, narg, arg)
{
  if (narg != iarg) error->all(FLERR,"Illegal compute inertia/molecule command");

  if (atom->molecular == 0)
    error->all(FLERR,"Compute inertia/molecule requires molecular atom style");

  array_flag = 1;
  size_array_cols = 6;
  extarray = 0;

  // setup molecule-based data

  nmolecules = molecules_in_group(idlo,idhi);
  size_array_rows = nmolecules;

  memory->create(massproc,nmolecules,"inertia/molecule:massproc");
  memory->create(masstotal,nmolecules,"inertia/molecule:masstotal");
  memory->create(com,nmolecules,3,"inertia/molecule:com");
  memory->create(comall,nmolecules,3,"inertia/molecule:comall");
  memory->create(inertia,nmolecules,6,"inertia/molecule:inertia");
  memory->create(inertiaall,nmolecules,6,"inertia/molecule:inertiaall");
  array = inertiaall;

  // compute masstotal for each molecule

  int *mask = atom->mask;
  int *molecule = atom->molecule;
  int *type = atom->type;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int nlocal = atom->nlocal;

  int i,imol;
  double massone;

  for (i = 0; i < nmolecules; i++) massproc[i] = 0.0;

  for (i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      if (rmass) massone = rmass[i];
      else massone = mass[type[i]];
      imol = molecule[i];
      if (molmap) imol = molmap[imol-idlo];
      else imol--;
      massproc[imol] += massone;
    }

  MPI_Allreduce(massproc,masstotal,nmolecules,MPI_DOUBLE,MPI_SUM,world);
}

/* ---------------------------------------------------------------------- */

ComputeInertiaMolecule::~ComputeInertiaMolecule()
{
  memory->destroy(massproc);
  memory->destroy(masstotal);
  memory->destroy(com);
  memory->destroy(comall);
  memory->destroy(inertia);
  memory->destroy(inertiaall);
}

/* ---------------------------------------------------------------------- */

void ComputeInertiaMolecule::init()
{
  int ntmp = molecules_in_group(idlo,idhi);
  if (ntmp != nmolecules)
    error->all(FLERR,"Molecule count changed in compute inertia/molecule");
}

/* ---------------------------------------------------------------------- */

void ComputeInertiaMolecule::compute_array()
{
  int i,j,imol;
  double dx,dy,dz,massone;
  double unwrap[3];

  invoked_array = update->ntimestep;

  double **x = atom->x;
  int *mask = atom->mask;
  int *molecule = atom->molecule;
  int *type = atom->type;
  tagint *image = atom->image;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int nlocal = atom->nlocal;

  // center-of-mass for each molecule

  for (i = 0; i < nmolecules; i++)
    com[i][0] = com[i][1] = com[i][2] = 0.0;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      if (rmass) massone = rmass[i];
      else massone = mass[type[i]];
      imol = molecule[i];
      if (molmap) imol = molmap[imol-idlo];
      else imol--;
      domain->unmap(x[i],image[i],unwrap);
      com[imol][0] += unwrap[0] * massone;
      com[imol][1] += unwrap[1] * massone;
      com[imol][2] += unwrap[2] * massone;
    }

  MPI_Allreduce(&com[0][0],&comall[0][0],3*nmolecules,
                MPI_DOUBLE,MPI_SUM,world);
  for (i = 0; i < nmolecules; i++) {
    comall[i][0] /= masstotal[i];
    comall[i][1] /= masstotal[i];
    comall[i][2] /= masstotal[i];
  }

  // inertia tensor for each molecule

  for (i = 0; i < nmolecules; i++)
    for (j = 0; j < 6; j++)
      inertia[i][j] = 0.0;

  for (i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      if (rmass) massone = rmass[i];
      else massone = mass[type[i]];
      imol = molecule[i];
      if (molmap) imol = molmap[imol-idlo];
      else imol--;
      domain->unmap(x[i],image[i],unwrap);
      dx = unwrap[0] - comall[imol][0];
      dy = unwrap[1] - comall[imol][1];
      dz = unwrap[2] - comall[imol][2];
      inertia[imol][0] += massone * (dy*dy + dz*dz);
      inertia[imol][1] += massone * (dx*dx + dz*dz);
      inertia[imol][2] += massone * (dx*dx + dy*dy);
      inertia[imol][3] -= massone * dx*dy;
      inertia[imol][4] -= massone * dy*dz;
      inertia[imol][5] -= massone * dx*dz;
    }

  MPI_Allreduce(&inertia[0][0],&inertiaall[0][0],6*nmolecules,
                MPI_DOUBLE,MPI_SUM,world);
}

/* ----------------------------------------------------------------------
   memory usage of local data
------------------------------------------------------------------------- */

double ComputeInertiaMolecule::memory_usage()
{
  double bytes = 2*nmolecules * sizeof(double);
  if (molmap) bytes += (idhi-idlo+1) * sizeof(int);
  bytes += 2*nmolecules*3 * sizeof(double);
  bytes += 2*nmolecules*6 * sizeof(double);
  return bytes;
}
