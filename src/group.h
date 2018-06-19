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

#ifndef LMP_GROUP_H
#define LMP_GROUP_H

#include <stdio.h>
#include "pointers.h"

namespace LAMMPS_NS {

class Group : protected Pointers {
 public:
  int ngroup;                  // # of defined groups
  char **names;                // name of each group
  int *bitmask;                // one-bit mask for each group
  int *inversemask;            // inverse mask for each group

  Group(class LAMMPS *);
  ~Group();
  void assign(int, char **);         // assign atoms to a group
  void create(const char *, int *);        // add flagged atoms to a group
  void init();                             // init function for general set up
  void set(const char *name, bool flag);   // set if all atoms belong to group or not
  int find(const char *);            // lookup name in list of groups
  void write_restart(FILE *);
  void read_restart(FILE *);

  bigint count(int);                       // count atoms in group
  bigint count(int,int);                   // count atoms in group & region
  bigint count_ms(int);                    // count atoms / multispheres in group
  bigint count_ms(int,int);                // count atoms / multispheres in group & region
  double mass(int);                        // total mass of atoms in group
  double mass(int,int);
  double charge(int);                      // total charge of atoms in group
  double charge(int,int);
  void bounds(int, double *);              // bounds of atoms in group
  void bounds(int, double *, int);
  void xcm(int, double, double *);         // center-of-mass coords of group
  void xcm(int, double, double *, int);
  void vcm(int, double, double *);         // center-of-mass velocity of group
  void vcm(int, double, double *, int);
  void fcm(int, double *);                 // total force on group
  void fcm(int, double *, int);
  double ke(int);                          // kinetic energy of group
  double ke(int, int);
  double gyration(int, double, double *);  // radius-of-gyration of group
  double gyration(int, double, double *, int);
  void angmom(int, double *, double *);    // angular momentum of group
  void angmom(int, double *, double *, int);
  void torque(int, double *, double *);    // torque on group
  void torque(int, double *, double *, int);
  void inertia(int, double *, double [3][3]);     // inertia tensor
  void inertia(int, double *, double [3][3], int);
  void omega(double *, double [3][3], double *);  // angular velocity

 private:
  int me;
  class FixMultisphere *fix_ms_; // holds multispheres, otherwise NULL

  int find_unused();
};

}

#endif

/* ERROR/WARNING messages:

E: Group command before simulation box is defined

The group command cannot be used before a read_data, read_restart, or
create_box command.

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Could not find group delete group ID

Self-explanatory.

E: Cannot delete group all

Self-explanatory.

E: Cannot delete group currently used by a fix

Self-explanatory.

E: Cannot delete group currently used by a compute

Self-explanatory.

E: Cannot delete group currently used by a dump

Self-explanatory.

E: Cannot delete group currently used by atom_modify first

Self-explanatory.

E: Too many groups

The maximum number of atom groups (including the "all" group) is
given by MAX_GROUP in group.cpp and is 32.

E: Group region ID does not exist

A region ID used in the group command does not exist.

E: Variable name for group does not exist

Self-explanatory.

E: Variable for group is invalid style

Only atom-style variables can be used.

E: Group ID does not exist

A group ID used in the group command does not exist.

*/
