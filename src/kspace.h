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

#ifndef LMP_KSPACE_H
#define LMP_KSPACE_H

#include "pointers.h"

#ifdef FFT_SINGLE
typedef float FFT_SCALAR;
#define MPI_FFT_SCALAR MPI_FLOAT
#else
typedef double FFT_SCALAR;
#define MPI_FFT_SCALAR MPI_DOUBLE
#endif

namespace LAMMPS_NS {

class KSpace : protected Pointers {
  friend class ThrOMP;
  friend class FixOMP;
 public:
  double energy;                 // accumulated energies
  double energy_1,energy_6;
  double virial[6];              // accumlated virial
  double *eatom,**vatom;         // accumulated per-atom energy/virial
  double e2group;                // accumulated group-group energy
  double f2group[3];             // accumulated group-group force
  int triclinic_support;         // 1 if supports triclinic geometries

  int ewaldflag;                 // 1 if a Ewald solver
  int pppmflag;                  // 1 if a PPPM solver
  int msmflag;                   // 1 if a MSM solver
  int dispersionflag;            // 1 if a LJ/dispersion solver
  int tip4pflag;                 // 1 if a TIP4P solver
  int dipoleflag;                // 1 if a dipole solver
  int differentiation_flag;
  int slabflag;
  double slab_volfactor;

  int order,order_6,order_allocated;
  double accuracy;                  // accuracy of KSpace solver (force units)
  double accuracy_absolute;         // user-specifed accuracy in force units
  double accuracy_relative;         // user-specified dimensionless accuracy
                                    // accurary = acc_rel * two_charge_force
  double two_charge_force;          // force in user units of two point
                                    // charges separated by 1 Angstrom

  double g_ewald,g_ewald_6;
  int nx_pppm,ny_pppm,nz_pppm;           // global FFT grid for Coulombics
  int nx_pppm_6,ny_pppm_6,nz_pppm_6;     // global FFT grid for dispersion
  int nx_msm_max,ny_msm_max,nz_msm_max;

  int group_group_enable;         // 1 if style supports group/group calculation

  unsigned int datamask;
  unsigned int datamask_ext;

  int compute_flag;               // 0 if skip compute()
  int fftbench;                   // 0 if skip FFT timing

  int stagger_flag;               // 1 if using staggered PPPM grids

  KSpace(class LAMMPS *, int, char **);
  virtual ~KSpace();
  void triclinic_check();
  void modify_params(int, char **);
  void *extract(const char *);
  void compute_dummy(int, int);

  // triclinic

  void x2lamdaT(double *, double *);
  void lamda2xT(double *, double *);
  void lamda2xvector(double *, double *);
  void kspacebbox(double, double *);

  // general child-class methods

  virtual void init() = 0;
  virtual void setup() = 0;
  virtual void setup_grid() {};
  virtual void compute(int, int) = 0;
  virtual void compute_group_group(int, int, int) {};

  virtual void pack_forward(int, FFT_SCALAR *, int, int *) {};
  virtual void unpack_forward(int, FFT_SCALAR *, int, int *) {};
  virtual void pack_reverse(int, FFT_SCALAR *, int, int *) {};
  virtual void unpack_reverse(int, FFT_SCALAR *, int, int *) {};

  virtual int timing(int, double &, double &) {return 0;}
  virtual int timing_1d(int, double &) {return 0;}
  virtual int timing_3d(int, double &) {return 0;}
  virtual double memory_usage() {return 0.0;}

/* ----------------------------------------------------------------------
   compute gamma for MSM and pair styles
   see Eq 4 from Parallel Computing 35 (2009) 164177
------------------------------------------------------------------------- */

  double gamma(const double &rho) const {
    if (rho <= 1.0) {
      const int split_order = order/2;
      const double rho2 = rho*rho;
      double g = gcons[split_order][0];
      double rho_n = rho2;
      for (int n=1; n<=split_order; n++) {
        g += gcons[split_order][n]*rho_n;
        rho_n *= rho2;
      }
      return g;
    } else
      return (1.0/rho);
  }

/* ----------------------------------------------------------------------
   compute the derivative of gamma for MSM and pair styles
   see Eq 4 from Parallel Computing 35 (2009) 164-177
------------------------------------------------------------------------- */

  double dgamma(const double &rho) const {
    if (rho <= 1.0) {
      const int split_order = order/2;
      const double rho2 = rho*rho;
      double dg = dgcons[split_order][0]*rho;
      double rho_n = rho*rho2;
      for (int n=1; n<split_order; n++) {
        dg += dgcons[split_order][n]*rho_n;
        rho_n *= rho2;
      }
      return dg;
    } else
      return (-1.0/rho/rho);
  }

  double **get_gcons() { return gcons; }
  double **get_dgcons() { return dgcons; }

 protected:
  int gridflag,gridflag_6;
  int gewaldflag,gewaldflag_6;
  int minorder,overlap_allowed;
  int adjust_cutoff_flag;
  int suffix_flag;                  // suffix compatibility flag
  double scale;
  double **gcons,**dgcons;          // accumulated per-atom energy/virial

  int evflag,evflag_atom;
  int eflag_either,eflag_global,eflag_atom;
  int vflag_either,vflag_global,vflag_atom;
  int maxeatom,maxvatom;

  int kewaldflag;                   // 1 if kspace range set for Ewald sum
  int kx_ewald,ky_ewald,kz_ewald;   // kspace settings for Ewald sum

  void pair_check();
  void ev_setup(int, int);
  double estimate_table_accuracy(double, double);
};

}

#endif

/* ERROR/WARNING messages:

E: KSpace style does not yet support triclinic geometries

UNDOCUMENTED

E: KSpace solver requires a pair style

No pair style is defined.

E: KSpace style is incompatible with Pair style

Setting a kspace style requires that a pair style with a long-range
Coulombic or dispersion component be used.

W: For better accuracy use 'pair_modify table 0'

The user-specified force accuracy cannot be achieved unless the table
feature is disabled by using 'pair_modify table 0'.

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Bad kspace_modify slab parameter

Kspace_modify value for the slab/volume keyword must be >= 2.0.

E: Bad kspace_modify kmax/ewald parameter

Kspace_modify values for the kmax/ewald keyword must be integers > 0

W: Kspace_modify slab param < 2.0 may cause unphysical behavior

The kspace_modify slab parameter should be larger to insure periodic
grids padded with empty space do not overlap.

*/
