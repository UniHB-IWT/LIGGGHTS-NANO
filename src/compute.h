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

#ifndef LMP_COMPUTE_H
#define LMP_COMPUTE_H

#include "pointers.h"

#define INVOKED_SCALAR  1<<0
#define INVOKED_VECTOR  1<<1
#define INVOKED_ARRAY   1<<2
#define INVOKED_PERATOM 1<<3
#define INVOKED_LOCAL   1<<4

namespace LAMMPS_NS {

class Compute : protected Pointers {
 public:
  char *id,*style;
  int igroup,groupbit;

  double scalar;            // computed global scalar
  double *vector;           // computed global vector
  double **array;           // computed global array
  double *vector_atom;      // computed per-atom vector
  double **array_atom;      // computed per-atom array
  double *vector_local;     // computed local vector
  double **array_local;     // computed local array

  int scalar_flag;          // 0/1 if compute_scalar() function exists
  int vector_flag;          // 0/1 if compute_vector() function exists
  int array_flag;           // 0/1 if compute_array() function exists
  int size_vector;          // length of global vector
  int size_array_rows;      // rows in global array
  int size_array_cols;      // columns in global array

  int peratom_flag;         // 0/1 if compute_peratom() function exists
  int size_peratom_cols;    // 0 = vector, N = columns in peratom array

  int local_flag;           // 0/1 if compute_local() function exists
  int size_local_rows;      // rows in local vector or array
  int size_local_cols;      // 0 = vector, N = columns in local array

  int extscalar;            // 0/1 if global scalar is intensive/extensive
  int extvector;            // 0/1/-1 if global vector is all int/ext/extlist
  int *extlist;             // list of 0/1 int/ext for each vec component
  int extarray;             // 0/1 if global array is all intensive/extensive

  int tempflag;       // 1 if Compute can be used as temperature
                      // must have both compute_scalar, compute_vector
  int pressflag;      // 1 if Compute can be used as pressure (uses virial)
                      // must have both compute_scalar, compute_vector
  int pressatomflag;  // 1 if Compute calculates per-atom virial
  int peflag;         // 1 if Compute calculates PE (uses Force energies)
  int peatomflag;     // 1 if Compute calculates per-atom PE

  int tempbias;       // 0/1 if Compute temp includes self/extra bias

  int timeflag;       // 1 if Compute stores list of timesteps it's called on
  int ntime;          // # of entries in time list
  int maxtime;        // max # of entries time list can hold
  bigint *tlist;      // list of timesteps the Compute is called on

  int invoked_flag;       // non-zero if invoked or accessed this step, 0 if not
  bigint invoked_scalar;  // last timestep on which compute_scalar() was invoked
  bigint invoked_vector;  // ditto for compute_vector()
  bigint invoked_array;   // ditto for compute_array()
  bigint invoked_peratom; // ditto for compute_peratom()
  bigint invoked_local;   // ditto for compute_local()

  double dof;         // degrees-of-freedom for temperature

  int comm_forward;   // size of forward communication (0 if none)
  int comm_reverse;   // size of reverse communication (0 if none)
  int dynamic_group_allow;  // 1 if can be used with dynamic group, else 0

  unsigned int datamask;
  unsigned int datamask_ext;
  int cudable;        // 1 if compute is CUDA-enabled

  Compute(class LAMMPS *lmp, int &iarg, int narg, char ** arg);
  virtual ~Compute();
  void modify_params(int, char **);
  void reset_extra_dof();

  virtual void post_create() {} 
  virtual void pre_delete(bool) {} 
  virtual void init() = 0;
  virtual void init_list(int, class NeighList *) {}
  virtual void setup() {}
  virtual double compute_scalar() {return 0.0;}
  virtual void compute_vector() {}
  virtual void compute_array() {}
  virtual void compute_peratom() {}
  virtual void compute_local() {}

  virtual int pack_comm(int, int *, double *, int, int *) {return 0;}
  virtual void unpack_comm(int, int, double *) {}
  virtual int pack_reverse_comm(int, int, double *) {return 0;}
  virtual void unpack_reverse_comm(int, int *, double *) {}

  virtual int dof_remove(int) {return 0;}
  virtual void remove_bias(int, double *) {}
  virtual void remove_bias_all() {}
  virtual void restore_bias(int, double *) {}
  virtual void restore_bias_all() {}

  virtual void reset_extra_compute_fix(const char *);

  void addstep(bigint);
  int matchstep(bigint);
  void clearstep();

  virtual double memory_usage() {return 0.0;}
  virtual int unsigned data_mask() {return datamask;}
  virtual int unsigned data_mask_ext() {return datamask_ext;}

  bool update_on_run_end()
  { return update_on_run_end_; }

 protected:
  int extra_dof;               // extra DOF for temperature computes
  int dynamic;                 // recount atoms for temperature computes
  int thermoflag;              // 1 if include fix PE for PE computes

  double vbias[3];             // stored velocity bias for one atom
  double **vbiasall;           // stored velocity bias for all atoms
  int maxbias;                 // size of vbiasall array

  int *molmap;                 // convert molecule ID to local index

  int molecules_in_group(int &, int &);

  inline int sbmask(int j) {
    return j >> SBBITS & 3;
  }

  // true if this compute is updated at the end of every run
  bool update_on_run_end_;
};

}

#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Compute ID must be alphanumeric or underscore characters

Self-explanatory.

E: Could not find compute group ID

Self-explanatory.

E: Compute does not allow an extra compute or fix to be reset

This is an internal LAMMPS error.  Please report it to the
developers.

W: Atom with molecule ID = 0 included in compute molecule group

The group used in a compute command that operates on moleclues
includes atoms with no molecule ID.  This is probably not what you
want.

W: One or more compute molecules has atoms not in group

The group used in a compute command that operates on moleclues does
not include all the atoms in some molecules.  This is probably not
what you want.

*/
