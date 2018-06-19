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

#ifdef ATOM_CLASS

AtomStyle(hybrid,AtomVecHybrid)

#else

#ifndef LMP_ATOM_VEC_HYBRID_H
#define LMP_ATOM_VEC_HYBRID_H

#include <stdio.h>
#include "atom_vec.h"

namespace LAMMPS_NS {

class AtomVecHybrid : public AtomVec {
 public:
  int nstyles;
  class AtomVec **styles;
  char **keywords;

  AtomVecHybrid(class LAMMPS *);
  ~AtomVecHybrid();
  void settings(int, char **);
  void init();
  void grow(int);
  void grow_reset();
  void copy(int, int, int);
  void clear_bonus();
  int pack_comm(int, int *, double *, int, int *);
  int pack_comm_vel(int, int *, double *, int, int *);
  void unpack_comm(int, int, double *);
  void unpack_comm_vel(int, int, double *);
  int pack_reverse(int, int, double *);
  void unpack_reverse(int, int *, double *);
  int pack_border(int, int *, double *, int, int *);
  int pack_border_vel(int, int *, double *, int, int *);
  void unpack_border(int, int, double *);
  void unpack_border_vel(int, int, double *);
  int pack_exchange(int, double *);
  int unpack_exchange(double *);
  int size_restart();
  int pack_restart(int, double *);
  int unpack_restart(double *);
  void write_restart_settings(FILE *);
  void read_restart_settings(FILE *);
  void create_atom(int, double *);
  void data_atom(double *, tagint, char **);
  int data_atom_hybrid(int, char **) {return 0;}
  void data_vel(int, char **);
  void pack_data(double **);
  void write_data(FILE *, int, double **);
  void pack_vel(double **);
  void write_vel(FILE *, int, double **);
  bigint memory_usage();

 private:
  int *tag,*type,*mask;
  tagint *image;
  double **x,**v,**f;
  double **omega,**angmom;

  int nallstyles;
  char **allstyles;

  void build_styles();
  int known_style(char *);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Atom style hybrid cannot have hybrid as an argument

Self-explanatory.

E: Atom style hybrid cannot use same atom style twice

Self-explanatory.

E: Per-processor system is too big

The number of owned atoms plus ghost atoms on a single
processor must fit in 32-bit integer.

E: Invalid atom ID in Atoms section of data file

Atom IDs must be positive integers.

E: Invalid atom type in Atoms section of data file

Atom types must range from 1 to specified # of types.

*/
