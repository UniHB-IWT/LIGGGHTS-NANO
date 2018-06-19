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
    This file is from LAMMPS, but has been modified. Copyright for
    modification:

    Copyright 2012-     DCS Computing GmbH, Linz
    Copyright 2009-2012 JKU Linz

    Copyright of original file:
    LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
    http://lammps.sandia.gov, Sandia National Laboratories
    Steve Plimpton, sjplimp@sandia.gov

    Copyright (2003) Sandia Corporation.  Under the terms of Contract
    DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
    certain rights in this software.  This software is distributed under
    the GNU General Public License.
------------------------------------------------------------------------- */
#ifdef FIX_CLASS

FixStyle(rigid,FixRigid)

#else

#ifndef LMP_FIX_RIGID_H
#define LMP_FIX_RIGID_H

#include "fix.h"

namespace LAMMPS_NS {

class FixRigid : public Fix {
 friend class FixTemplateMultisphere;
 friend class FixGravity;
 friend class PairGran;
 friend class FixWallGran;

 public:
  FixRigid(class LAMMPS *, int, char **);
  virtual ~FixRigid();
  virtual int setmask();
  virtual void init();
  virtual void setup(int);
  virtual void initial_integrate(int);
  void post_force(int);
  virtual void final_integrate();
  void initial_integrate_respa(int, int, int);
  void final_integrate_respa(int, int);
  void write_restart_file(char *);
  virtual double compute_scalar();
  virtual int modify_param(int, char **) {return 0;}

  double memory_usage();
  void grow_arrays(int);
  void copy_arrays(int, int, int);
  void set_arrays(int);
  int pack_exchange(int, double *);
  int unpack_exchange(int, double *);

  void pre_neighbor();
  int dof(int);
  void deform(int);
  void reset_dt();
  void zero_momentum();
  void zero_rotation();
  virtual void *extract(const char*, int &);
  double extract_ke();
  double extract_erotational();
  double compute_array(int, int);

 protected:
  int me,nprocs;
  double dtv,dtf,dtq;
  double *step_respa;
  int triclinic;
  double MINUSPI,TWOPI;

  char *infile;             // file to read rigid body attributes from
  int rstyle;               // SINGLE,MOLECULE,GROUP
  int firstflag;            // 1 for first-time setup of rigid bodies

  int dimension;            // # of dimensions
  int nbody;                // # of rigid bodies
  int *nrigid;              // # of atoms in each rigid body
  int *mol2body;            // convert mol-ID to rigid body index
  int *body2mol;            // convert rigid body index to mol-ID
  int maxmol;               // size of mol2body = max mol-ID

  int *body;                // which body each atom is part of (-1 if none)
  double **displace;        // displacement of each atom in body coords
  double *masstotal;        // total mass of each rigid body
  double **xcm;             // coords of center-of-mass of each rigid body
  double **vcm;             // velocity of center-of-mass of each
  double **fcm;             // force on center-of-mass of each
  double **inertia;         // 3 principal components of inertia of each
  double **ex_space,**ey_space,**ez_space;
                            // principal axes of each in space coords
  double **angmom;          // angular momentum of each in space coords
  double **omega;           // angular velocity of each in space coords
  double **torque;          // torque on each rigid body in space coords
  double **quat;            // quaternion of each rigid body
  tagint *imagebody;        // image flags of xcm of each rigid body
  double **fflag;           // flag for on/off of center-of-mass force
  double **tflag;           // flag for on/off of center-of-mass torque
  double **langextra;       // Langevin thermostat forces and torques

  double **sum,**all;       // work vectors for each rigid body
  int **remapflag;          // PBC remap flags for each rigid body

  int extended;             // 1 if any particles have extended attributes
  int orientflag;           // 1 if particles store spatial orientation
  int dorientflag;          // 1 if particles store dipole orientation

  int *eflags;              // flags for extended particles
  double **orient;          // orientation vector of particle wrt rigid body
  double **dorient;         // orientation of dipole mu wrt rigid body

  double tfactor;           // scale factor on temperature of rigid bodies
  int langflag;             // 0/1 = no/yes Langevin thermostat

  int tstat_flag;           // NVT settings
  double t_start,t_stop,t_target;
  double t_period,t_freq;
  int t_chain,t_iter,t_order;

  int pstat_flag;           // NPT settings
  double p_start[3],p_stop[3];
  double p_period[3],p_freq[3];
  int p_flag[3];
  int pcouple,pstyle;
  int p_chain;
  int allremap;              // remap all atoms
  int dilate_group_bit;      // mask for dilation group
  char *id_dilate;           // group name to dilate

  class RanMars *random;
  class AtomVecEllipsoid *avec_ellipsoid;
  class AtomVecLine *avec_line;

  int POINT,SPHERE,ELLIPSOID,LINE,TRIANGLE,DIPOLE;   // bitmasks for eflags
  int OMEGA,ANGMOM,TORQUE;

  void no_squish_rotate(int, double *, double *, double *, double) const;
  void set_xv();
  void set_v();
  void setup_bodies_static();
  void setup_bodies_dynamic();
  void readfile(int, double *, double **, int *);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Fix rigid molecule requires atom attribute molecule

Self-explanatory.

E: Could not find fix rigid group ID

A group ID used in the fix rigid command does not exist.

E: One or more atoms belong to multiple rigid bodies

Two or more rigid bodies defined by the fix rigid command cannot
contain the same atom.

E: No rigid bodies defined

The fix specification did not end up defining any rigid bodies.

E: Fix rigid z force cannot be on for 2d simulation

Self-explanatory.

E: Fix rigid xy torque cannot be on for 2d simulation

Self-explanatory.

E: Fix rigid langevin period must be > 0.0

Self-explanatory.

E: Fix rigid npt/nph dilate group ID does not exist

Self-explanatory.
E: One or zero atoms in rigid body

Any rigid body defined by the fix rigid command must contain 2 or more
atoms.

W: More than one fix rigid

It is not efficient to use fix rigid more than once.

E: Rigid fix must come before NPT/NPH fix

NPT/NPH fix must be defined in input script after all rigid fixes,
else the rigid fix contribution to the pressure virial is
incorrect.
W: Cannot count rigid body degrees-of-freedom before bodies are fully initialized

UNDOCUMENTED

W: Computing temperature of portions of rigid bodies

The group defined by the temperature compute does not encompass all
the atoms in one or more rigid bodies, so the change in
degrees-of-freedom for the atoms in those partial rigid bodies will
not be accounted for.

E: Fix rigid atom has non-zero image flag in a non-periodic dimension

Image flags for non-periodic dimensions should not be set.

E: Insufficient Jacobi rotations for rigid body

Eigensolve for rigid body was not sufficiently accurate.

E: Fix rigid: Bad principal moments

The principal moments of inertia computed for a rigid body
are not within the required tolerances.

E: Cannot open fix rigid infile %s

The specified file cannot be opened.  Check that the path and name are
correct.

E: Unexpected end of fix rigid file

A read operation from the file failed.

E: Incorrect rigid body format in fix rigid file

The number of fields per line is not what expected.

E: Invalid rigid body ID in fix rigid file

The ID does not match the number or an existing ID of rigid bodies
that are defined by the fix rigid command.

*/
