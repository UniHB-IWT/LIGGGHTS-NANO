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
    (if not contributing author is listed, this file has been contributed
    by the core developer)

    Christoph Kloss (DCS Computing GmbH, JKU Linz)
    Arno Mayrhofer (DCS Computing GmbH)

    Copyright 2012-     DCS Computing GmbH, Linz
    Copyright 2009-2012 JKU Linz
------------------------------------------------------------------------- */

#ifdef MESHMODULE_CLASS

MeshModuleStyle(stress,MeshModuleStress)

#else

#ifndef LMP_MESH_MODULE_STRESS_H
#define LMP_MESH_MODULE_STRESS_H

#include "mesh_module.h"

namespace LAMMPS_NS
{
  class MeshModuleStress : public MeshModule
  {
      public:

        MeshModuleStress(LAMMPS *lmp, int &iarg_, int narg, char **arg, FixMeshSurface *fix_mesh);
        virtual ~MeshModuleStress();

        virtual void post_create_pre_restart();
        virtual void post_create();

        virtual void init();
        virtual void setup(int vflag) {}
        virtual int setmask();

        virtual void pre_force(int vflag);
        virtual void final_integrate();

        virtual double compute_vector(int n);

        virtual void add_particle_contribution(int ip, double *frc,
                            double *delta, int iTri, double *v_wall);

        void add_global_external_contribution(double *frc);
        void add_global_external_contribution(double *frc,double *trq);

        // inline access

        inline bool trackWear()
        { return wear_flag_; }

        inline bool trackStress()
        { return stress_flag_; }

        inline void f_total(double *_f)
        { vectorCopy3D(f_total_,_f); }

        inline double f_total(int i)
        { return f_total_[i]; }

        inline double *f_total_ptr()
        { return f_total_; }

        inline double f_total_mag()
        { return vectorMag3D(f_total_); }

        inline void torque_total(double *_t)
        { vectorCopy3D(torque_total_,_t); }

        inline double torque_total(int i)
        { return torque_total_[i]; }

        inline double *torque_total_ptr()
        { return torque_total_; }

        inline int get_num_vector_components() const
        { return 9; }

        inline void set_p_ref(double *_p_ref)
        { p_ref_.set(0,_p_ref); }

        inline double get_wear_increment(int i)
        { return (*wear_increment_)(i); }

        void set_store_wear_increment(const bool set)
        { store_wear_increment_ = set; }

      protected:

        // STRESS
        // total force and total torque
        double f_total_[3], torque_total_[3]; 

        // inline access

        inline double* f(int i)
        { return (*f_)(i); }

        inline double& sigma_n(int i)
        { return (*sigma_n_)(i); }

        inline double& sigma_t(int i)
        { return (*sigma_t_)(i); }

        inline double& wear(int i)
        { return (*wear_)(i); }

        inline double& wear_step(int i)
        { return (*wear_step_)(i); }

        inline double& wear_increment(int i)
        { return (*wear_increment_)(i); }

        inline double p_ref(int i)
        { return p_ref_(0)[i]; }

      private:

        // inititalization fcts
        void regStress();
        void regWear();
        void zeroizeStress();
        void zeroizeWear();

        void calc_total_force();
        void add_gravity();

        // STRESS
        //
        bool stress_flag_;

        // buffer variables for correct output
        double f_total_old_[3], torque_total_old_[3];
        bool updatedStresses_;

        // stress flag in FixMeshSurface
        // reference point, total force and total torque
        VectorContainer<double,3> &p_ref_;

        // per-element force and torque
        VectorContainer<double,3> *f_;
        ScalarContainer<double> *sigma_n_;
        ScalarContainer<double> *sigma_t_;

        // WEAR

        // flag for wear model and Finnie constant
        int wear_flag_;
        double const* const* k_finnie_;
        ScalarContainer<double> *wear_;
        ScalarContainer<double> *wear_step_;
        ScalarContainer<double> *wear_increment_;
        bool store_wear_increment_;
  };

} /* namespace LAMMPS_NS */

#endif
#endif
