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

#ifdef COMMAND_CLASS

CommandStyle(displace_atoms,DisplaceAtoms)

#else

#ifndef LMP_DISPLACE_ATOMS_H
#define LMP_DISPLACE_ATOMS_H

#include "pointers.h"

namespace LAMMPS_NS {

class DisplaceAtoms : protected Pointers {
 public:
  DisplaceAtoms(class LAMMPS *);
  void command(int, char **);

 private:
  int scaleflag;
  void options(int, char **);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Displace_atoms command before simulation box is defined

The displace_atoms command cannot be used before a read_data,
read_restart, or create_box command.

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Cannot displace_atoms after reading restart file with per-atom info

This is because the restart file info cannot be migrated with the
atoms.  You can get around this by performing a 0-timestep run which
will assign the restart file info to actual atoms.

E: Could not find displace_atoms group ID

Group ID used in the displace_atoms command does not exist.

E: Invalid displace_atoms rotate axis for 2d

Axis must be in z direction.

E: Zero length rotation vector with displace_atoms

Self-explanatory.

W: Lost atoms via displace_atoms: original %ld current %ld

The command options you have used caused atoms to be lost.

U: Use of displace_atoms with undefined lattice

Must use lattice command with displace_atoms command if units option
is set to lattice.

*/
