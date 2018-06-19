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

#ifdef FIX_CLASS

FixStyle(move,FixMove)

#else

#ifndef LMP_FIX_MOVE_H
#define LMP_FIX_MOVE_H

#include <stdio.h>
#include "fix.h"

namespace LAMMPS_NS {

class FixMove : public Fix {
 public:
  FixMove(class LAMMPS *, int, char **);
  ~FixMove();
  int setmask();
  void init();
  void initial_integrate(int);
  void final_integrate();
  void initial_integrate_respa(int, int, int);
  void final_integrate_respa(int, int);

  double memory_usage();
  void write_restart(FILE *);
  void restart(char *);
  void grow_arrays(int);
  void copy_arrays(int, int, int);
  void set_arrays(int);
  int pack_exchange(int, double *);
  int unpack_exchange(int, double *);
  int pack_restart(int, double *);
  void unpack_restart(int, int);
  int maxsize_restart();
  int size_restart(int);

  void reset_dt();

 private:
  char *xvarstr,*yvarstr,*zvarstr,*vxvarstr,*vyvarstr,*vzvarstr;
  int mstyle;
  int vxflag,vyflag,vzflag,axflag,ayflag,azflag;
  double vx,vy,vz,ax,ay,az;
  double period,omega_rotate;
  double point[3],axis[3],runit[3];
  double dt,dtv,dtf;
  int xvar,yvar,zvar,vxvar,vyvar,vzvar;
  int xvarstyle,yvarstyle,zvarstyle,vxvarstyle,vyvarstyle,vzvarstyle;
  int omega_flag,nlevels_respa;
  int time_origin;

  double **xoriginal;         // original coords of atoms
  int displaceflag,velocityflag;
  int maxatom;
  double **displace,**velocity;
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Fix move cannot set linear z motion for 2d problem

Self-explanatory.

E: Fix move cannot set wiggle z motion for 2d problem

Self-explanatory.

E: Fix move cannot rotate aroung non z-axis for 2d problem

Self-explanatory.

E: Fix move cannot define z or vz variable for 2d problem

Self-explanatory.

W: Fix move does not update angular momentum

Atoms store this quantity, but fix move does not (yet) update it.

W: Fix move does not update quaternions

Atoms store this quantity, but fix move does not (yet) update it.

E: Zero length rotation vector with fix move

Self-explanatory.

E: Variable name for fix move does not exist

Self-explanatory.

E: Variable for fix move is invalid style

Only equal-style variables can be used.

E: Cannot add atoms to fix move variable

Atoms can not be added afterwards to this fix option.

E: Resetting timestep is not allowed with fix move

This is because fix move is moving atoms based on elapsed time.

U: Use of fix move with undefined lattice

Must use lattice command with fix move command if units option is
set to lattice.

*/
