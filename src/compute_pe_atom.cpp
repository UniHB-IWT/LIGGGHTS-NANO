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
#include "compute_pe_atom.h"
#include "atom.h"
#include "update.h"
#include "comm.h"
#include "force.h"
#include "pair.h"
#include "bond.h"
#include "angle.h"
#include "dihedral.h"
#include "improper.h"
#include "kspace.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputePEAtom::ComputePEAtom(LAMMPS *lmp, int &iarg, int narg, char **arg) :
  Compute(lmp, iarg, narg, arg)
{
  if (narg < iarg) error->all(FLERR,"Illegal compute pe/atom command");

  peratom_flag = 1;
  size_peratom_cols = 0;
  peatomflag = 1;
  timeflag = 1;
  comm_reverse = 1;

  if (narg == iarg) {
    pairflag = 1;
    bondflag = angleflag = dihedralflag = improperflag = 1;
    kspaceflag = 1;
  } else {
    pairflag = 0;
    bondflag = angleflag = dihedralflag = improperflag = 0;
    kspaceflag = 0;
    while (iarg < narg) {
      if (strcmp(arg[iarg],"pair") == 0) pairflag = 1;
      else if (strcmp(arg[iarg],"bond") == 0) bondflag = 1;
      else if (strcmp(arg[iarg],"angle") == 0) angleflag = 1;
      else if (strcmp(arg[iarg],"dihedral") == 0) dihedralflag = 1;
      else if (strcmp(arg[iarg],"improper") == 0) improperflag = 1;
      else if (strcmp(arg[iarg],"kspace") == 0) kspaceflag = 1;
      else error->all(FLERR,"Illegal compute pe/atom command");
      iarg++;
    }
  }

  nmax = 0;
  energy = NULL;
}

/* ---------------------------------------------------------------------- */

ComputePEAtom::~ComputePEAtom()
{
  memory->destroy(energy);
}

/* ---------------------------------------------------------------------- */

void ComputePEAtom::compute_peratom()
{
  int i;

  invoked_peratom = update->ntimestep;
  if (update->eflag_atom != invoked_peratom)
    error->all(FLERR,"Per-atom energy was not tallied on needed timestep");

  // grow local energy array if necessary
  // needs to be atom->nmax in length

  if (atom->nmax > nmax) {
    memory->destroy(energy);
    nmax = atom->nmax;
    memory->create(energy,nmax,"pe/atom:energy");
    vector_atom = energy;
  }

  // npair includes ghosts if either newton flag is set
  //   b/c some bonds/dihedrals call pair::ev_tally with pairwise info
  // nbond includes ghosts if newton_bond is set
  // ntotal includes ghosts if either newton flag is set
  // KSpace includes ghosts if tip4pflag is set

  int nlocal = atom->nlocal;
  int npair = nlocal;
  int nbond = nlocal;
  int ntotal = nlocal;
  int nkspace = nlocal;
  if (force->newton) npair += atom->nghost;
  if (force->newton_bond) nbond += atom->nghost;
  if (force->newton) ntotal += atom->nghost;
  if (force->kspace && force->kspace->tip4pflag) nkspace += atom->nghost;

  // clear local energy array

  for (i = 0; i < ntotal; i++) energy[i] = 0.0;

  // add in per-atom contributions from each force

  if (pairflag && force->pair) {
    double *eatom = force->pair->eatom;
    for (i = 0; i < npair; i++) energy[i] += eatom[i];
  }

  if (bondflag && force->bond) {
    double *eatom = force->bond->eatom;
    for (i = 0; i < nbond; i++) energy[i] += eatom[i];
  }

  if (angleflag && force->angle) {
    double *eatom = force->angle->eatom;
    for (i = 0; i < nbond; i++) energy[i] += eatom[i];
  }

  if (dihedralflag && force->dihedral) {
    double *eatom = force->dihedral->eatom;
    for (i = 0; i < nbond; i++) energy[i] += eatom[i];
  }

  if (improperflag && force->improper) {
    double *eatom = force->improper->eatom;
    for (i = 0; i < nbond; i++) energy[i] += eatom[i];
  }

  if (kspaceflag && force->kspace) {
    double *eatom = force->kspace->eatom;
    for (i = 0; i < nkspace; i++) energy[i] += eatom[i];
  }

  // communicate ghost energy between neighbor procs

  if (force->newton || (force->kspace && force->kspace->tip4pflag))
    comm->reverse_comm_compute(this);

  // zero energy of atoms not in group
  // only do this after comm since ghost contributions must be included

  int *mask = atom->mask;

  for (i = 0; i < nlocal; i++)
    if (!(mask[i] & groupbit)) energy[i] = 0.0;
}

/* ---------------------------------------------------------------------- */

int ComputePEAtom::pack_reverse_comm(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) buf[m++] = energy[i];
  return 1;
}

/* ---------------------------------------------------------------------- */

void ComputePEAtom::unpack_reverse_comm(int n, int *list, double *buf)
{
  int i,j,m;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    energy[j] += buf[m++];
  }
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double ComputePEAtom::memory_usage()
{
  double bytes = nmax * sizeof(double);
  return bytes;
}
