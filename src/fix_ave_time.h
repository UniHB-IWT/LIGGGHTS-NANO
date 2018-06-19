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

FixStyle(ave/time,FixAveTime)

#else

#ifndef LMP_FIX_AVE_TIME_H
#define LMP_FIX_AVE_TIME_H

#include <stdio.h>
#include "fix.h"

namespace LAMMPS_NS {

class FixAveTime : public Fix {
 public:
  FixAveTime(class LAMMPS *, int, char **);
  ~FixAveTime();
  int setmask();
  void init();
  void setup(int);
  void end_of_step();
  double compute_scalar();
  double compute_vector(int);
  double compute_array(int,int);
  void reset_timestep(bigint);

 private:
  int me,nvalues;
  int nrepeat,nfreq,irepeat;
  bigint nvalid;
  int *which,*argindex,*value2index,*offcol;
  char **ids;
  FILE *fp;
  int nrows;

  int ave,nwindow,nsum,startstep,mode;
  int noff,overwrite;
  int *offlist;
  char *title1,*title2,*title3;
  long filepos;

  int norm,iwindow,window_limit;
  double *vector;
  double *vector_total;
  double **vector_list;
  double *column;
  double **array;
  double **array_total;
  double ***array_list;

  void invoke_scalar(bigint);
  void invoke_vector(bigint);
  void options(int, char **);
  void allocate_values(int);
  bigint nextvalid();
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Compute ID for fix ave/time does not exist

Self-explanatory.

E: Fix ID for fix ave/time does not exist

Self-explanatory.

E: Invalid fix ave/time off column

Self-explantory.

E: Fix ave/time compute does not calculate a scalar

Self-explantory.

E: Fix ave/time compute does not calculate a vector

Self-explantory.

E: Fix ave/time compute vector is accessed out-of-range

The index for the vector is out of bounds.

E: Fix ave/time compute does not calculate an array

Self-explanatory.

E: Fix ave/time compute array is accessed out-of-range

An index for the array is out of bounds.

E: Fix ave/time fix does not calculate a scalar

Self-explanatory.

E: Fix ave/time fix does not calculate a vector

Self-explanatory.

E: Fix ave/time fix vector is accessed out-of-range

The index for the vector is out of bounds.

E: Fix for fix ave/time not computed at compatible time

Fixes generate their values on specific timesteps.  Fix ave/time
is requesting a value on a non-allowed timestep.

E: Fix ave/time fix does not calculate an array

Self-explanatory.

E: Fix ave/time fix array is accessed out-of-range

An index for the array is out of bounds.

E: Variable name for fix ave/time does not exist

Self-explanatory.

E: Fix ave/time variable is not equal-style variable

Self-explanatory.

E: Fix ave/time cannot use variable with vector mode

Variables produce scalar values.

E: Fix ave/time columns are inconsistent lengths

Self-explanatory.

E: Fix ave/time cannot set output array intensive/extensive from these inputs

One of more of the vector inputs has individual elements which are
flagged as intensive or extensive.  Such an input cannot be flagged as
all intensive/extensive when turned into an array by fix ave/time.

E: Cannot open fix ave/time file %s

The specified file cannot be opened.  Check that the path and name are
correct.

E: Fix ave/time missed timestep

You cannot reset the timestep to a value beyond where the fix
expects to next perform averaging.

*/
