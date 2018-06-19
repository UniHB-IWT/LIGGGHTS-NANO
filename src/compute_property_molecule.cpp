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
#include "compute_property_molecule.h"
#include "atom.h"
#include "update.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputePropertyMolecule::ComputePropertyMolecule(LAMMPS *lmp, int &iarg, int narg, char **arg) :
  Compute(lmp, iarg, narg, arg)
{
  if (narg < iarg+1) error->all(FLERR,"Illegal compute property/molecule command");

  if (atom->molecular == 0)
    error->all(FLERR,"Compute property/molecule requires molecular atom style");

  nvalues = narg - iarg;

  pack_choice = new FnPtrPack[nvalues];

  int i;
  const int arg_offset = iarg;
  for (; iarg < narg; iarg++) {
    i = iarg-arg_offset;

    if (strcmp(arg[iarg],"mol") == 0)
      pack_choice[i] = &ComputePropertyMolecule::pack_mol;
    else if (strcmp(arg[iarg],"count") == 0)
      pack_choice[i] = &ComputePropertyMolecule::pack_count;
    else error->all(FLERR,"Invalid keyword in compute property/molecule command");
  }

  // setup molecule-based data

  nmolecules = molecules_in_group(idlo,idhi);

  vector = NULL;
  array = NULL;

  if (nvalues == 1) {
    memory->create(vector,nmolecules,"property/molecule:vector");
    vector_flag = 1;
    size_vector = nmolecules;
    extvector = 0;
  } else {
    memory->create(array,nmolecules,nvalues,"property/molecule:array");
    array_flag = 1;
    size_array_rows = nmolecules;
    size_array_cols = nvalues;
    extarray = 0;
  }

  // fill vector or array with molecule values

  if (nvalues == 1) {
    buf = vector;
    (this->*pack_choice[0])(0);
  } else {
    if (array) buf = &array[0][0];
    for (int n = 0; n < nvalues; n++)
      (this->*pack_choice[n])(n);
  }
}

/* ---------------------------------------------------------------------- */

ComputePropertyMolecule::~ComputePropertyMolecule()
{
  delete [] pack_choice;
  memory->destroy(vector);
  memory->destroy(array);
}

/* ---------------------------------------------------------------------- */

void ComputePropertyMolecule::init()
{
  int ntmp = molecules_in_group(idlo,idhi);
  if (ntmp != nmolecules)
    error->all(FLERR,"Molecule count changed in compute property/molecule");
}

/* ---------------------------------------------------------------------- */

void ComputePropertyMolecule::compute_vector()
{
  invoked_vector = update->ntimestep;
}

/* ---------------------------------------------------------------------- */

void ComputePropertyMolecule::compute_array()
{
  invoked_array = update->ntimestep;
}

/* ----------------------------------------------------------------------
   memory usage of local data
------------------------------------------------------------------------- */

double ComputePropertyMolecule::memory_usage()
{
  double bytes = nmolecules*nvalues * sizeof(double);
  if (molmap) bytes += (idhi-idlo+1) * sizeof(int);
  return bytes;
}

/* ----------------------------------------------------------------------
   one method for every keyword compute property/molecule can output
   the atom property is packed into buf starting at n with stride nvalues
   customize a new keyword by adding a method
------------------------------------------------------------------------- */

void ComputePropertyMolecule::pack_mol(int n)
{
  for (int m = idlo; m <= idhi; m++)
    if (molmap == NULL || molmap[m-idlo] >= 0) {
      buf[n] = m;
      n += nvalues;
    }
}

/* ---------------------------------------------------------------------- */

void ComputePropertyMolecule::pack_count(int n)
{
  int i,m,imol;

  int *count_one = new int[nmolecules];
  for (m = 0; m < nmolecules; m++) count_one[m] = 0;

  int *molecule = atom->molecule;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      imol = molecule[i];
      if (molmap) imol = molmap[imol-idlo];
      else imol--;
      count_one[imol]++;
    }

  int *count_all = new int[nmolecules];
  MPI_Allreduce(count_one,count_all,nmolecules,MPI_INT,MPI_SUM,world);

  for (m = 0; m < nmolecules; m++)
    if (molmap == NULL || molmap[m] >= 0) {
      buf[n] = count_all[m];
      n += nvalues;
    }

  delete [] count_one;
  delete [] count_all;
}
