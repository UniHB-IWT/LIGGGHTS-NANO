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

/* ----------------------------------------------------------------------
   Contributing authors:
   Richard Berger (JKU Linz)
------------------------------------------------------------------------- */

#include <stdio.h>
#include <string.h>
#include "modify.h"
#include "style_compute.h"
#include "style_fix.h"
#include "atom.h"
#include "comm.h"
#include "fix.h"
#include "compute.h"
#include "group.h"
#include "update.h"
#include "domain.h"
#include "memory.h"
#include "error.h"
#include "force.h"
#include <map>
#include <string>

using namespace LAMMPS_NS;
using namespace FixConst;

#define DELTA 4

#define BIG 1.0e20
#define NEXCEPT 4       // change when add to exceptions in add_fix()

/* ---------------------------------------------------------------------- */

Modify::Modify(LAMMPS *lmp) : Pointers(lmp)
{
  nfix = maxfix = 0;
  n_pre_initial_integrate = n_initial_integrate = n_post_integrate = 0;
  n_pre_exchange = n_pre_neighbor = 0;
  n_pre_force = n_post_force = 0;
  n_iterate_implicitly = n_pre_final_integrate = 0; 
  n_final_integrate = n_end_of_step = n_thermo_energy = 0;
  n_initial_integrate_respa = n_post_integrate_respa = 0;
  n_pre_force_respa = n_post_force_respa = n_final_integrate_respa = 0;
  n_min_pre_exchange = n_min_pre_force = n_min_post_force = n_min_energy = 0;
  n_min_pre_neighbor = 0;

  fix = NULL;
  fmask = NULL;
  list_pre_initial_integrate = list_initial_integrate = list_post_integrate = NULL;
  list_pre_exchange = list_pre_neighbor = NULL;
  list_pre_force = list_post_force = NULL;
  list_iterate_implicitly = list_pre_final_integrate = NULL; 
  list_final_integrate = list_end_of_step = NULL;
  list_thermo_energy = NULL;
  list_initial_integrate_respa = list_post_integrate_respa = NULL;
  list_pre_force_respa = list_post_force_respa = NULL;
  list_final_integrate_respa = NULL;
  list_min_pre_exchange = list_min_pre_neighbor = NULL;
  list_min_pre_force = list_min_post_force = NULL;
  list_min_energy = NULL;

  end_of_step_every = NULL;

  list_timeflag = NULL;

  nfix_restart_global = 0;
  id_restart_global = style_restart_global = state_restart_global = NULL;
  nfix_restart_peratom = 0;
  id_restart_peratom = style_restart_peratom = NULL;
  index_restart_peratom = NULL;

  ncompute = maxcompute = 0;
  compute = NULL;

  timing = 0;

  // fill map with fixes listed in style_fix.h

  fix_map = new std::map<std::string,FixCreator>();

#define FIX_CLASS
#define FixStyle(key,Class) \
  (*fix_map)[#key] = &fix_creator<Class>;
#include "style_fix.h"
#undef FixStyle
#undef FIX_CLASS

  // fill map with computes listed in style_compute.h

  compute_map = new std::map<std::string,ComputeCreator>();

#define COMPUTE_CLASS
#define ComputeStyle(key,Class) \
  (*compute_map)[#key] = &compute_creator<Class>;
#include "style_compute.h"
#undef ComputeStyle
#undef COMPUTE_CLASS
}

/* ---------------------------------------------------------------------- */

Modify::~Modify()
{
  // delete all fixes
  // do it via delete_fix() so callbacks in Atom are also updated correctly

  while (nfix) delete_fix(fix[0]->id);
  memory->sfree(fix);
  memory->destroy(fmask);

  // delete all computes

  for (int i = ncompute-1; i >= 0; i--) delete compute[i];
  memory->sfree(compute);

  delete [] list_pre_initial_integrate;
  delete [] list_initial_integrate;
  delete [] list_post_integrate;
  delete [] list_pre_exchange;
  delete [] list_pre_neighbor;
  delete [] list_pre_force;
  delete [] list_post_force;
  delete [] list_pre_final_integrate;
  delete [] list_final_integrate;
  delete [] list_iterate_implicitly; 
  delete [] list_end_of_step;
  delete [] list_thermo_energy;
  delete [] list_initial_integrate_respa;
  delete [] list_post_integrate_respa;
  delete [] list_pre_force_respa;
  delete [] list_post_force_respa;
  delete [] list_final_integrate_respa;
  delete [] list_min_pre_exchange;
  delete [] list_min_pre_neighbor;
  delete [] list_min_pre_force;
  delete [] list_min_post_force;
  delete [] list_min_energy;

  delete [] end_of_step_every;
  delete [] list_timeflag;

  restart_deallocate();
  delete compute_map;
  delete fix_map;
}

/* ----------------------------------------------------------------------
   initialize all fixes and computes
------------------------------------------------------------------------- */

void Modify::init()
{
  int i,j;

  // delete storage of restart info since it is not valid after 1st run

  restart_deallocate();

  // create lists of fixes to call at each stage of run

  list_init(PRE_INITIAL_INTEGRATE,n_pre_initial_integrate,list_pre_initial_integrate);
  list_init(INITIAL_INTEGRATE,n_initial_integrate,list_initial_integrate);
  list_init(POST_INTEGRATE,n_post_integrate,list_post_integrate);
  list_init_pre_exchange(PRE_EXCHANGE,n_pre_exchange,list_pre_exchange);
  //list_init(PRE_EXCHANGE,n_pre_exchange,list_pre_exchange);
  list_init(PRE_NEIGHBOR,n_pre_neighbor,list_pre_neighbor);
  list_init(PRE_FORCE,n_pre_force,list_pre_force);
  list_init(POST_FORCE,n_post_force,list_post_force);
  list_init(PRE_FINAL_INTEGRATE,n_pre_final_integrate,list_pre_final_integrate);
  list_init(FINAL_INTEGRATE,n_final_integrate,list_final_integrate);
  list_init(ITERATE_IMPLICITLY,n_iterate_implicitly,list_iterate_implicitly); 
  list_init_end_of_step(END_OF_STEP,n_end_of_step,list_end_of_step);
  list_init_thermo_energy(THERMO_ENERGY,n_thermo_energy,list_thermo_energy);

  list_init(INITIAL_INTEGRATE_RESPA,
            n_initial_integrate_respa,list_initial_integrate_respa);
  list_init(POST_INTEGRATE_RESPA,
            n_post_integrate_respa,list_post_integrate_respa);
  list_init(POST_FORCE_RESPA,
            n_post_force_respa,list_post_force_respa);
  list_init(PRE_FORCE_RESPA,
            n_pre_force_respa,list_pre_force_respa);
  list_init(FINAL_INTEGRATE_RESPA,
            n_final_integrate_respa,list_final_integrate_respa);

  list_init(MIN_PRE_EXCHANGE,n_min_pre_exchange,list_min_pre_exchange);
  list_init(MIN_PRE_FORCE,n_min_pre_force,list_min_pre_force);
  list_init(MIN_POST_FORCE,n_min_post_force,list_min_post_force);
  list_init(MIN_ENERGY,n_min_energy,list_min_energy);

  // init each fix
  // not sure if now needs to come before compute init
  // used to b/c temperature computes called fix->dof() in their init,
  // and fix rigid required its own init before its dof() could be called,
  // but computes now do their DOF in setup()

  for (i = 0; i < nfix; i++) fix[i]->init();

  // set global flag if any fix has its restart_pbc flag set

  restart_pbc_any = 0;
  for (i = 0; i < nfix; i++)
    if (fix[i]->restart_pbc) restart_pbc_any = 1;

  // create list of computes that store invocation times

  list_init_compute();

  // init each compute
  // set invoked_scalar,vector,etc to -1 to force new run to re-compute them
  // add initial timestep to all computes that store invocation times
  //   since any of them may be invoked by initial thermo
  // do not clear out invocation times stored within a compute,
  //   b/c some may be holdovers from previous run, like for ave fixes

  for (i = 0; i < ncompute; i++) {
    compute[i]->init();
    compute[i]->invoked_scalar = -1;
    compute[i]->invoked_vector = -1;
    compute[i]->invoked_array = -1;
    compute[i]->invoked_peratom = -1;
    compute[i]->invoked_local = -1;
  }
  addstep_compute_all(update->ntimestep);

  // warn if any particle is time integrated more than once

  int nlocal = atom->nlocal;
  int *mask = atom->mask;

  int *flag = new int[nlocal];
  for (i = 0; i < nlocal; i++) flag[i] = 0;

  int groupbit;
  for (i = 0; i < nfix; i++) {
    if (fix[i]->time_integrate == 0) continue;
    groupbit = fix[i]->groupbit;
    for (j = 0; j < nlocal; j++)
      if (mask[j] & groupbit) flag[j]++;
  }

  int check = 0;
  for (i = 0; i < nlocal; i++)
    if (flag[i] > 1) check = 1;

  delete [] flag;

  int checkall;
  MPI_Allreduce(&check,&checkall,1,MPI_INT,MPI_SUM,world);
  if (comm->me == 0 && checkall)
    error->warning(FLERR,
                   "One or more atoms are time integrated more than once");
}

/* ----------------------------------------------------------------------
   setup for run, calls setup() of all fixes and computes
   called from Verlet, RESPA, Min
------------------------------------------------------------------------- */

void Modify::setup(int vflag)
{
  
  int nlocal = atom->nlocal;
  int *mask = atom->mask;

  for (int i = 0; i < nfix; i++)
  {
     if (!fix[i]->recent_restart && fix[i]->just_created && fix[i]->create_attribute)
     {
         fix[i]->just_created = 0;
         fix[i]->pre_set_arrays();
         for(int j = 0; j < nlocal; j++)
           if(mask[j] & fix[i]->groupbit) fix[i]->set_arrays(j);
     }
     else if(fix[i]->just_created) fix[i]->just_created = 0;
     fix[i]->recent_restart = 0;
  }

  // compute setup needs to come before fix setup
  // b/c NH fixes need use DOF of temperature computes

  for (int i = 0; i < ncompute; i++) compute[i]->setup();

  if (update->whichflag == 1)
    call_method_on_fixes(&Fix::setup, vflag);
  else if (update->whichflag == 2)
    call_method_on_fixes(&Fix::min_setup, vflag);
}

/* ----------------------------------------------------------------------
   setup pre_exchange call, only for fixes that define pre_exchange
   called from Verlet, RESPA, Min, and WriteRestart with whichflag = 0
------------------------------------------------------------------------- */

void Modify::setup_pre_exchange()
{
  if (update->whichflag <= 1)
    call_method_on_fixes(&Fix::setup_pre_exchange, list_pre_exchange, n_pre_exchange);
  else if (update->whichflag == 2)
    call_method_on_fixes(&Fix::min_setup_pre_exchange, list_min_pre_exchange, n_min_pre_exchange);
}

/* ----------------------------------------------------------------------
   setup pre_neighbor call, only for fixes that define pre_neighbor
   called from Verlet, RESPA
------------------------------------------------------------------------- */

void Modify::setup_pre_neighbor()
{
  if (update->whichflag == 1)
    call_method_on_fixes(&Fix::setup_pre_neighbor, list_pre_neighbor, n_pre_neighbor);
  else if (update->whichflag == 2)
    call_method_on_fixes(&Fix::min_setup_pre_neighbor, list_min_pre_neighbor, n_min_pre_neighbor);
}

/* ----------------------------------------------------------------------
   setup pre_force call, only for fixes that define pre_force
   called from Verlet, RESPA, Min
------------------------------------------------------------------------- */

void Modify::setup_pre_force(int vflag)
{
  if (update->whichflag == 1) {
    call_method_on_fixes(&Fix::setup_pre_force, vflag, list_pre_force, n_pre_force);
  } else if (update->whichflag == 2) {
    call_method_on_fixes(&Fix::min_setup_pre_force, vflag, list_min_pre_force, n_min_pre_force);
  }
}

/* ----------------------------------------------------------------------
   pre_intial_integrate call, only for relevant fixes
------------------------------------------------------------------------- */

void Modify::pre_initial_integrate()
{
  call_method_on_fixes(&Fix::pre_initial_integrate, list_pre_initial_integrate, n_pre_initial_integrate);
}

/* ----------------------------------------------------------------------
   1st half of integrate call, only for relevant fixes
------------------------------------------------------------------------- */

void Modify::initial_integrate(int vflag)
{
  
  call_method_on_fixes(&Fix::initial_integrate, vflag, list_initial_integrate, n_initial_integrate);
}

/* ----------------------------------------------------------------------
   post_integrate call, only for relevant fixes
------------------------------------------------------------------------- */

void Modify::post_integrate()
{
  call_method_on_fixes(&Fix::post_integrate, list_post_integrate, n_post_integrate);
}

/* ----------------------------------------------------------------------
   pre_exchange call, only for relevant fixes
------------------------------------------------------------------------- */

void Modify::pre_exchange()
{
  
  call_method_on_fixes(&Fix::pre_exchange, list_pre_exchange, n_pre_exchange);
}

/* ----------------------------------------------------------------------
   pre_neighbor call, only for relevant fixes
------------------------------------------------------------------------- */

void Modify::pre_neighbor()
{
  
  call_method_on_fixes(&Fix::pre_neighbor, list_pre_neighbor, n_pre_neighbor);
}

/* ----------------------------------------------------------------------
   pre_force call, only for relevant fixes
------------------------------------------------------------------------- */

void Modify::pre_force(int vflag)
{
  
  call_method_on_fixes(&Fix::pre_force, vflag, list_pre_force, n_pre_force);
}

/* ----------------------------------------------------------------------
   post_force call, only for relevant fixes
------------------------------------------------------------------------- */

void Modify::post_force(int vflag)
{
  
  call_method_on_fixes(&Fix::post_force, vflag, list_post_force, n_post_force);
}

/* ----------------------------------------------------------------------
   check convergence, only for relevant implicit integration
------------------------------------------------------------------------- */

bool Modify::iterate_implicitly()
{
  for (int i = 0; i < n_iterate_implicitly; i++)
    if (fix[list_iterate_implicitly[i]]->iterate_implicitly())
        return true;

  return false;
}

/* ----------------------------------------------------------------------
   pre_final_integrate call, only for relevant fixes
------------------------------------------------------------------------- */

void Modify::pre_final_integrate()
{
  call_method_on_fixes(&Fix::pre_final_integrate, list_pre_final_integrate, n_pre_final_integrate);
}

/* ----------------------------------------------------------------------
   2nd half of integrate call, only for relevant fixes
------------------------------------------------------------------------- */

void Modify::final_integrate()
{
  call_method_on_fixes(&Fix::final_integrate, list_final_integrate, n_final_integrate);
}

/* ----------------------------------------------------------------------
   end-of-timestep call, only for relevant fixes
   only call fix->end_of_step() on timesteps that are multiples of nevery
------------------------------------------------------------------------- */

void Modify::end_of_step()
{
  if(timing) {
    for (int i = 0; i < n_end_of_step; i++) {
      if (update->ntimestep % end_of_step_every[i] == 0) {
        const int ifix = list_end_of_step[i];
        fix[ifix]->begin_time_recording();
        fix[ifix]->end_of_step();
        fix[ifix]->end_time_recording();
      }
    }
  }
  else
  {
    for (int i = 0; i < n_end_of_step; i++) {
      if (update->ntimestep % end_of_step_every[i] == 0) {
        fix[list_end_of_step[i]]->end_of_step();
      }
    }
  }
}

/* ----------------------------------------------------------------------
   thermo energy call, only for relevant fixes
   called by Thermo class
   compute_scalar() is fix call to return energy
------------------------------------------------------------------------- */

double Modify::thermo_energy()
{
  double energy = 0.0;
  if(timing) {
    for (int i = 0; i < n_thermo_energy; i++) {
      const int ifix = list_thermo_energy[i];
      fix[ifix]->begin_time_recording();
      energy += fix[ifix]->compute_scalar();
      fix[ifix]->end_time_recording();
    }
  }
  else
  {
    for (int i = 0; i < n_thermo_energy; i++) {
      energy += fix[list_thermo_energy[i]]->compute_scalar();
    }
  }
  return energy;
}

/* ----------------------------------------------------------------------
   post_run call
------------------------------------------------------------------------- */

void Modify::post_run()
{
  call_method_on_fixes(&Fix::post_run);
}

/* ----------------------------------------------------------------------
   setup rRESPA pre_force call, only for relevant fixes
------------------------------------------------------------------------- */

void Modify::setup_pre_force_respa(int vflag, int ilevel)
{
  call_respa_method_on_fixes(&Fix::setup_pre_force_respa, vflag, ilevel,
      list_pre_force_respa, n_pre_force_respa);
}

/* ----------------------------------------------------------------------
   1st half of rRESPA integrate call, only for relevant fixes
------------------------------------------------------------------------- */

void Modify::initial_integrate_respa(int vflag, int ilevel, int iloop)
{
  call_respa_method_on_fixes(&Fix::initial_integrate_respa, vflag, ilevel, iloop,
      list_initial_integrate_respa, n_initial_integrate_respa);
}

/* ----------------------------------------------------------------------
   rRESPA post_integrate call, only for relevant fixes
------------------------------------------------------------------------- */

void Modify::post_integrate_respa(int ilevel, int iloop)
{
  call_respa_method_on_fixes(&Fix::post_integrate_respa, ilevel, iloop,
      list_post_integrate_respa, n_post_integrate_respa);
}

/* ----------------------------------------------------------------------
   rRESPA pre_force call, only for relevant fixes
------------------------------------------------------------------------- */

void Modify::pre_force_respa(int vflag, int ilevel, int iloop)
{
  call_respa_method_on_fixes(&Fix::pre_force_respa, vflag, ilevel, iloop,
      list_pre_force_respa, n_pre_force_respa);
}

/* ----------------------------------------------------------------------
   rRESPA post_force call, only for relevant fixes
------------------------------------------------------------------------- */

void Modify::post_force_respa(int vflag, int ilevel, int iloop)
{
  call_respa_method_on_fixes(&Fix::post_force_respa, vflag, ilevel, iloop,
      list_post_force_respa, n_post_force_respa);
}

/* ----------------------------------------------------------------------
   2nd half of rRESPA integrate call, only for relevant fixes
------------------------------------------------------------------------- */

void Modify::final_integrate_respa(int ilevel, int iloop)
{
  call_respa_method_on_fixes(&Fix::final_integrate_respa, ilevel, iloop,
      list_final_integrate_respa, n_final_integrate_respa);
}

/* ----------------------------------------------------------------------
   minimizer pre-exchange call, only for relevant fixes
------------------------------------------------------------------------- */

void Modify::min_pre_exchange()
{
  call_method_on_fixes(&Fix::min_pre_exchange, list_min_pre_exchange, n_min_pre_exchange);
}

/* ----------------------------------------------------------------------
   minimizer pre-neighbor call, only for relevant fixes
------------------------------------------------------------------------- */

void Modify::min_pre_neighbor()
{
  for (int i = 0; i < n_min_pre_neighbor; i++)
    fix[list_min_pre_neighbor[i]]->min_pre_neighbor();
}

/* ----------------------------------------------------------------------
   minimizer pre-force call, only for relevant fixes
------------------------------------------------------------------------- */

void Modify::min_pre_force(int vflag)
{
  call_method_on_fixes(&Fix::min_pre_force, vflag, list_min_pre_force, n_min_pre_force);
}

/* ----------------------------------------------------------------------
   minimizer force adjustment call, only for relevant fixes
------------------------------------------------------------------------- */

void Modify::min_post_force(int vflag)
{
  call_method_on_fixes(&Fix::min_post_force, vflag, list_min_post_force, n_min_post_force);
}

/* ----------------------------------------------------------------------
   minimizer energy/force evaluation, only for relevant fixes
   return energy and forces on extra degrees of freedom
------------------------------------------------------------------------- */

double Modify::min_energy(double *fextra)
{
  int ifix,index;

  index = 0;
  double eng = 0.0;
  for (int i = 0; i < n_min_energy; i++) {
    ifix = list_min_energy[i];
    eng += fix[ifix]->min_energy(&fextra[index]);
    index += fix[ifix]->min_dof();
  }
  return eng;
}

/* ----------------------------------------------------------------------
   store current state of extra dof, only for relevant fixes
------------------------------------------------------------------------- */

void Modify::min_store()
{
  call_method_on_fixes(&Fix::min_store, list_min_energy, n_min_energy);
}

/* ----------------------------------------------------------------------
   mange state of extra dof on a stack, only for relevant fixes
------------------------------------------------------------------------- */

void Modify::min_clearstore()
{
  call_method_on_fixes(&Fix::min_clearstore, list_min_energy, n_min_energy);
}

void Modify::min_pushstore()
{
  call_method_on_fixes(&Fix::min_pushstore, list_min_energy, n_min_energy);
}

void Modify::min_popstore()
{
  call_method_on_fixes(&Fix::min_popstore, list_min_energy, n_min_energy);
}

/* ----------------------------------------------------------------------
   displace extra dof along vector hextra, only for relevant fixes
------------------------------------------------------------------------- */

void Modify::min_step(double alpha, double *hextra)
{
  int ifix,index;

  index = 0;
  for (int i = 0; i < n_min_energy; i++) {
    ifix = list_min_energy[i];
    fix[ifix]->min_step(alpha,&hextra[index]);
    index += fix[ifix]->min_dof();
  }
}

/* ----------------------------------------------------------------------
   compute max allowed step size along vector hextra, only for relevant fixes
------------------------------------------------------------------------- */

double Modify::max_alpha(double *hextra)
{
  int ifix,index;

  double alpha = BIG;
  index = 0;
  for (int i = 0; i < n_min_energy; i++) {
    ifix = list_min_energy[i];
    double alpha_one = fix[ifix]->max_alpha(&hextra[index]);
    alpha = MIN(alpha,alpha_one);
    index += fix[ifix]->min_dof();
  }
  return alpha;
}

/* ----------------------------------------------------------------------
   extract extra dof for minimization, only for relevant fixes
------------------------------------------------------------------------- */

int Modify::min_dof()
{
  int ndof = 0;
  for (int i = 0; i < n_min_energy; i++)
    ndof += fix[list_min_energy[i]]->min_dof();
  return ndof;
}

/* ----------------------------------------------------------------------
   reset reference state of fix, only for relevant fixes
------------------------------------------------------------------------- */

int Modify::min_reset_ref()
{
  int itmp,itmpall;
  itmpall = 0;
  for (int i = 0; i < n_min_energy; i++) {
    itmp = fix[list_min_energy[i]]->min_reset_ref();
    if (itmp) itmpall = 1;
  }
  return itmpall;
}

/* ----------------------------------------------------------------------
   add a new fix or replace one with same ID
------------------------------------------------------------------------- */

void Modify::add_fix(int narg, char **arg, char *suffix)
{
  
  if(update->ntimestep_reset_since_last_run && fix_restart_in_progress())
    error->all(FLERR,"In case of restart, command 'reset_timestep' must come immediately before 'run'");

  if (narg < 3) error->all(FLERR,"Illegal fix command");

  // cannot define fix before box exists unless style is in exception list
  // don't like this way of checking for exceptions by adding fixes to list,
  //   but can't think of better way
  // too late if instantiate fix, then check flag set in fix constructor,
  //   since some fixes access domain settings in their constructor
  // change NEXCEPT above when add new fix to this list

  const char *exceptions[NEXCEPT] = {"GPU","OMP","property/atom","cmap"};

  if (domain->box_exist == 0) {
    int m;
    for (m = 0; m < NEXCEPT; m++)
      if (strcmp(arg[2],exceptions[m]) == 0) break;
    if (m == NEXCEPT)
      error->all(FLERR,"Fix command before simulation box is defined");
  }

  // check group ID

  int igroup = group->find(arg[1]);
  if (igroup == -1) error->all(FLERR,"Could not find fix group ID");

  // if fix ID exists:
  //   set newflag = 0 so create new fix in same location in fix list
  //   error if new style does not match old style
  //     since can't replace it (all when-to-invoke ptrs would be invalid)
  //   warn if new group != old group
  //   delete old fix, but do not call update_callback(),
  //     since will replace this fix and thus other fix locs will not change
  //   set ptr to NULL in case new fix scans list of fixes,
  //     e.g. scan will occur in add_callback() if called by new fix
  // if fix ID does not exist:
  //   set newflag = 1 so create new fix
  //   extend fix and fmask lists as necessary

  int ifix,newflag;
  for (ifix = 0; ifix < nfix; ifix++)
    if (strcmp(arg[0],fix[ifix]->id) == 0) break;

  if (ifix < nfix) {
    newflag = 0;
    
    if (strncmp(fix[ifix]->style,"insert/",7) == 0)
      error->all(FLERR,"Using a fix insert/* ID twice, which is not possible. Please use different ID");
    if (strcmp(arg[2],fix[ifix]->style) != 0)
      error->all(FLERR,"Replacing a fix, but new style != old style");
    if (fix[ifix]->igroup != igroup && comm->me == 0)
      error->warning(FLERR,"Replacing a fix, but new group != old group");
    
    fix[ifix]->pre_delete(true);
    delete fix[ifix];
    fix[ifix] = NULL;
  } else {
    newflag = 1;
    if (nfix == maxfix) {
      maxfix += DELTA;
      fix = (Fix **) memory->srealloc(fix,maxfix*sizeof(Fix *),"modify:fix");
      memory->grow(fmask,maxfix,"modify:fmask");
    }
  }

  // create the Fix
  // try first with suffix appended

  fix[ifix] = NULL;

  if (suffix && lmp->suffix_enable) {
    char estyle[256];
    sprintf(estyle,"%s/%s",arg[2],suffix);

    if (fix_map->find(estyle) != fix_map->end()) {
      FixCreator fix_creator = (*fix_map)[estyle];
      fix[ifix] = fix_creator(lmp,narg,arg);
    }
  }

  if (fix[ifix] == NULL && fix_map->find(arg[2]) != fix_map->end()) {
    FixCreator fix_creator = (*fix_map)[arg[2]];
    fix[ifix] = fix_creator(lmp,narg,arg);
  }

  if (fix[ifix] == NULL && strncmp(arg[2], "mesh/surface", 12) == 0)
  {
    FixCreator fix_creator = (*fix_map)["mesh/surface"];
    fix[ifix] = fix_creator(lmp, narg, arg);
  }

  if (fix[ifix] == NULL){ 
    char * errmsg = new char[30+strlen(arg[2])]; 
    sprintf(errmsg,"Invalid fix style: \"%s\"",arg[2]);
    error->all(FLERR,errmsg);
    delete [] errmsg; 
  } 

  // set fix mask values and increment nfix (if new)

  fmask[ifix] = fix[ifix]->setmask();
  if (newflag) nfix++;

  fix[ifix]->post_create_pre_restart(); 

  // check if Fix is in restart_global list
  // if yes, pass state info to the Fix so it can reset itself

  for (int i = 0; i < nfix_restart_global; i++)
  {
    
    if (strcmp(id_restart_global[i],fix[ifix]->id) == 0 &&
          (strcmp(style_restart_global[i],fix[ifix]->style) == 0 ||
            (fix[ifix]->accepts_restart_data_from_style && strcmp(style_restart_global[i],fix[ifix]->accepts_restart_data_from_style) == 0)
          )
       )
      {
          fix[ifix]->restart(state_restart_global[i]);
          fix[ifix]->recent_restart = 1; 
          if (comm->me == 0) {
            char *str = (char *) ("Resetting global state of Fix %s Style %s "
                                  "from restart file info\n");
            if (screen) fprintf(screen,str,fix[ifix]->id,fix[ifix]->style);
            if (logfile) fprintf(logfile,str,fix[ifix]->id,fix[ifix]->style);
      }
    }
  }

  // check if Fix is in restart_peratom list
  // if yes, loop over atoms so they can extract info from atom->extra array

  for (int i = 0; i < nfix_restart_peratom; i++)
  {
    if (strcmp(id_restart_peratom[i],fix[ifix]->id) == 0 &&
          (strcmp(style_restart_peratom[i],fix[ifix]->style) == 0 ||
            (fix[ifix]->accepts_restart_data_from_style && strcmp(style_restart_peratom[i],fix[ifix]->accepts_restart_data_from_style) == 0)
          )
       )
    {
      for (int j = 0; j < atom->nlocal; j++)
        fix[ifix]->unpack_restart(j,index_restart_peratom[i]);
      fix[ifix]->recent_restart = 1; 
      fix[ifix]->restart_reset = 1;
      if (comm->me == 0) {
        char *str = (char *) ("Resetting per-atom state of Fix %s Style %s "
                     "from restart file info\n");
        if (screen) fprintf(screen,str,fix[ifix]->id,fix[ifix]->style);
        if (logfile) fprintf(logfile,str,fix[ifix]->id,fix[ifix]->style);

      }
    }
  }

  fix[ifix]->post_create(); 
}

/* ----------------------------------------------------------------------
   one instance per fix in style_fix.h
------------------------------------------------------------------------- */

template <typename T>
Fix *Modify::fix_creator(LAMMPS *lmp, int narg, char **arg)
{
  return new T(lmp,narg,arg);
}

/* ----------------------------------------------------------------------
   modify a Fix's parameters
------------------------------------------------------------------------- */

void Modify::modify_fix(int narg, char **arg)
{
  if (narg < 2) error->all(FLERR,"Illegal fix_modify command");

  // lookup Fix ID

  int ifix;
  for (ifix = 0; ifix < nfix; ifix++)
    if (strcmp(arg[0],fix[ifix]->id) == 0) break;
  if (ifix == nfix) error->all(FLERR,"Could not find fix_modify ID");

  fix[ifix]->modify_params(narg-1,&arg[1]);
}

/* ----------------------------------------------------------------------
   delete a Fix from list of Fixes
   Atom class must update indices in its list of callbacks to fixes
------------------------------------------------------------------------- */

void Modify::delete_fix(const char *id, bool unfixflag) 
{
  int ifix = find_fix(id);
  if (ifix < 0) error->all(FLERR,"Could not find fix ID to delete");

  fix[ifix]->pre_delete(unfixflag);

  delete fix[ifix];
  fix[ifix] = NULL; 
  atom->update_callback(ifix);

  // move other Fixes and fmask down in list one slot
  
  for (int i = ifix+1; i < nfix; i++) fix[i-1] = fix[i];
  for (int i = ifix+1; i < nfix; i++) fmask[i-1] = fmask[i];
  nfix--;
}

/* ----------------------------------------------------------------------
   find a fix by ID
   return index of fix or -1 if not found
------------------------------------------------------------------------- */

int Modify::find_fix(const char *id)
{
  int ifix;
  for (ifix = 0; ifix < nfix; ifix++)
    if (strcmp(id,fix[ifix]->id) == 0) break;
  if (ifix == nfix) return -1;
  return ifix;
}

/* ----------------------------------------------------------------------
   add a new compute
------------------------------------------------------------------------- */

void Modify::add_compute(int narg, char **arg, char *suffix)
{
  if (narg < 3) error->all(FLERR,"Illegal compute command");

  // error check

  for (int icompute = 0; icompute < ncompute; icompute++)
    if (strcmp(arg[0],compute[icompute]->id) == 0)
      error->all(FLERR,"Reuse of compute ID");

  // extend Compute list if necessary

  if (ncompute == maxcompute) {
    maxcompute += DELTA;
    compute = (Compute **)
      memory->srealloc(compute,maxcompute*sizeof(Compute *),"modify:compute");
  }

  // create the Compute
  // try first with suffix appended

  compute[ncompute] = NULL;

  if (suffix && lmp->suffix_enable)
  {
      char estyle[256];
      sprintf(estyle,"%s/%s",arg[2],suffix);
      if (compute_map->find(estyle) != compute_map->end())
      {
          ComputeCreator compute_creator = (*compute_map)[estyle];
          int iarg = 0;
          compute[ncompute] = compute_creator(lmp, iarg, narg, arg);
      }
  }

  if (compute[ncompute] == NULL &&
      compute_map->find(arg[2]) != compute_map->end())
  {
      int iarg = 0;
      ComputeCreator compute_creator = (*compute_map)[arg[2]];
      compute[ncompute] = compute_creator(lmp, iarg, narg, arg);
  }

  if (compute[ncompute] == NULL) error->all(FLERR,"Invalid compute style");

  ncompute++;

  compute[ncompute-1]->post_create(); 
}

/* ----------------------------------------------------------------------
   one instance per compute in style_compute.h
------------------------------------------------------------------------- */

template <typename T>
Compute *Modify::compute_creator(LAMMPS *lmp, int iarg, int narg, char **arg)
{
  return new T(lmp, iarg, narg,arg);
}

/* ----------------------------------------------------------------------
   modify a Compute's parameters
------------------------------------------------------------------------- */

void Modify::modify_compute(int narg, char **arg)
{
  if (narg < 2) error->all(FLERR,"Illegal compute_modify command");

  // lookup Compute ID

  int icompute;
  for (icompute = 0; icompute < ncompute; icompute++)
    if (strcmp(arg[0],compute[icompute]->id) == 0) break;
  if (icompute == ncompute)
    error->all(FLERR,"Could not find compute_modify ID");

  compute[icompute]->modify_params(narg-1,&arg[1]);
}

/* ----------------------------------------------------------------------
   delete a Compute from list of Computes
------------------------------------------------------------------------- */

void Modify::delete_compute(const char *id,bool uncomputeflag)
{
  int icompute = find_compute(id);
  if (icompute < 0) error->all(FLERR,"Could not find compute ID to delete");

  compute[icompute]->pre_delete(uncomputeflag);

  delete compute[icompute];

  // move other Computes down in list one slot

  for (int i = icompute+1; i < ncompute; i++) compute[i-1] = compute[i];
  ncompute--;
}

/* ----------------------------------------------------------------------
   find a compute by ID
   return index of compute or -1 if not found
------------------------------------------------------------------------- */

int Modify::find_compute(const char *id)
{
  int icompute;
  for (icompute = 0; icompute < ncompute; icompute++)
    if (strcmp(id,compute[icompute]->id) == 0) break;
  if (icompute == ncompute) return -1;
  return icompute;
}

/* ----------------------------------------------------------------------
   clear invoked flag of all computes
   called everywhere that computes are used, before computes are invoked
   invoked flag used to avoid re-invoking same compute multiple times
   and to flag computes that store invocation times as having been invoked
------------------------------------------------------------------------- */

void Modify::clearstep_compute()
{
  for (int icompute = 0; icompute < ncompute; icompute++)
    compute[icompute]->invoked_flag = 0;
}

/* ----------------------------------------------------------------------
   loop over computes that store invocation times
   if its invoked flag set on this timestep, schedule next invocation
   called everywhere that computes are used, after computes are invoked
------------------------------------------------------------------------- */

void Modify::addstep_compute(bigint newstep)
{
  for (int icompute = 0; icompute < n_timeflag; icompute++)
    if (compute[list_timeflag[icompute]]->invoked_flag)
      compute[list_timeflag[icompute]]->addstep(newstep);
}

/* ----------------------------------------------------------------------
   loop over all computes
   schedule next invocation for those that store invocation times
   called when not sure what computes will be needed on newstep
   do not loop only over n_timeflag, since may not be set yet
------------------------------------------------------------------------- */

void Modify::addstep_compute_all(bigint newstep)
{
  for (int icompute = 0; icompute < ncompute; icompute++)
    if (compute[icompute]->timeflag) compute[icompute]->addstep(newstep);
}

/* ----------------------------------------------------------------------
   write to restart file for all Fixes with restart info
   (1) fixes that have global state
   (2) fixes that store per-atom quantities
------------------------------------------------------------------------- */

void Modify::write_restart(FILE *fp)
{
  int me = comm->me;

  int count = 0;
  for (int i = 0; i < nfix; i++)
    if (fix[i]->restart_global) count++;

  if (me == 0) fwrite(&count,sizeof(int),1,fp);

  int n;
  for (int i = 0; i < nfix; i++)
    if (fix[i]->restart_global) {
      if (me == 0) {
        
        n = strlen(fix[i]->id) + 1;
        fwrite(&n,sizeof(int),1,fp);
        fwrite(fix[i]->id,sizeof(char),n,fp);
        n = strlen(fix[i]->style) + 1;
        fwrite(&n,sizeof(int),1,fp);
        fwrite(fix[i]->style,sizeof(char),n,fp);
      }
      fix[i]->write_restart(fp);
    }

  count = 0;
  for (int i = 0; i < nfix; i++)
    if (fix[i]->restart_peratom) count++;

  if (me == 0) fwrite(&count,sizeof(int),1,fp);

  for (int i = 0; i < nfix; i++)
    if (fix[i]->restart_peratom) {
      int maxsize_restart = fix[i]->maxsize_restart();
      if (me == 0) {
        
        n = strlen(fix[i]->id) + 1;
        fwrite(&n,sizeof(int),1,fp);
        fwrite(fix[i]->id,sizeof(char),n,fp);
        n = strlen(fix[i]->style) + 1;
        fwrite(&n,sizeof(int),1,fp);
        fwrite(fix[i]->style,sizeof(char),n,fp);
        fwrite(&maxsize_restart,sizeof(int),1,fp);
      }
    }
}

/* ----------------------------------------------------------------------
   read in restart file data on all previously defined Fixes with restart info
   (1) fixes that have global state
   (2) fixes that store per-atom quantities
   return maxsize of extra info that will be stored with any atom
------------------------------------------------------------------------- */

int Modify::read_restart(FILE *fp)
{
  // nfix_restart_global = # of restart entries with global state info

  int me = comm->me;
  if (me == 0) fread(&nfix_restart_global,sizeof(int),1,fp);
  MPI_Bcast(&nfix_restart_global,1,MPI_INT,0,world);

  // allocate space for each entry

  if (nfix_restart_global) {
    id_restart_global = new char*[nfix_restart_global];
    style_restart_global = new char*[nfix_restart_global];
    state_restart_global = new char*[nfix_restart_global];
  }

  // read each entry and Bcast to all procs
  // each entry has id string, style string, chunk of state data

  int n;
  for (int i = 0; i < nfix_restart_global; i++) {
    if (me == 0) fread(&n,sizeof(int),1,fp);
    MPI_Bcast(&n,1,MPI_INT,0,world);
    id_restart_global[i] = new char[n];
    if (me == 0) fread(id_restart_global[i],sizeof(char),n,fp);
    MPI_Bcast(id_restart_global[i],n,MPI_CHAR,0,world);

    if (me == 0) fread(&n,sizeof(int),1,fp);
    MPI_Bcast(&n,1,MPI_INT,0,world);
    style_restart_global[i] = new char[n];
    if (me == 0) fread(style_restart_global[i],sizeof(char),n,fp);
    MPI_Bcast(style_restart_global[i],n,MPI_CHAR,0,world);

    if (me == 0) fread(&n,sizeof(int),1,fp);
    MPI_Bcast(&n,1,MPI_INT,0,world);
    state_restart_global[i] = new char[n];
    if (me == 0) fread(state_restart_global[i],sizeof(char),n,fp);
    MPI_Bcast(state_restart_global[i],n,MPI_CHAR,0,world);
  }

  // nfix_restart_peratom = # of restart entries with peratom info

  int maxsize = 0;

  if (me == 0) fread(&nfix_restart_peratom,sizeof(int),1,fp);
  MPI_Bcast(&nfix_restart_peratom,1,MPI_INT,0,world);

  // allocate space for each entry

  if (nfix_restart_peratom) {
    id_restart_peratom = new char*[nfix_restart_peratom];
    style_restart_peratom = new char*[nfix_restart_peratom];
    index_restart_peratom = new int[nfix_restart_peratom];
  }

  // read each entry and Bcast to all procs
  // each entry has id string, style string, maxsize of one atom's data
  // set index = which set of extra data this fix represents

  for (int i = 0; i < nfix_restart_peratom; i++) {
    if (me == 0) fread(&n,sizeof(int),1,fp);
    MPI_Bcast(&n,1,MPI_INT,0,world);
    id_restart_peratom[i] = new char[n];
    if (me == 0) fread(id_restart_peratom[i],sizeof(char),n,fp);
    MPI_Bcast(id_restart_peratom[i],n,MPI_CHAR,0,world);

    if (me == 0) fread(&n,sizeof(int),1,fp);
    MPI_Bcast(&n,1,MPI_INT,0,world);
    style_restart_peratom[i] = new char[n];
    if (me == 0) fread(style_restart_peratom[i],sizeof(char),n,fp);
    MPI_Bcast(style_restart_peratom[i],n,MPI_CHAR,0,world);

    if (me == 0) fread(&n,sizeof(int),1,fp);
    MPI_Bcast(&n,1,MPI_INT,0,world);
    maxsize += n;

    index_restart_peratom[i] = i;

  }

  return maxsize;
}

/* ----------------------------------------------------------------------
   delete all lists of restart file Fix info
------------------------------------------------------------------------- */

void Modify::restart_deallocate()
{
  
  int n_ms = n_fixes_style("multisphere");
  bool have_ms_in_restart = false;

  if (nfix_restart_global) {
    for (int i = 0; i < nfix_restart_global; i++) {
      if(strncmp(style_restart_global[i],"multisphere",11) == 0)
        have_ms_in_restart = true;
      delete [] id_restart_global[i];
      delete [] style_restart_global[i];
      delete [] state_restart_global[i];
    }
    delete [] id_restart_global;
    delete [] style_restart_global;
    delete [] state_restart_global;
  }

  if (nfix_restart_peratom) {
    for (int i = 0; i < nfix_restart_peratom; i++) {
      if(strncmp(style_restart_peratom[i],"multisphere",11) == 0)
        have_ms_in_restart = true;
      delete [] id_restart_peratom[i];
      delete [] style_restart_peratom[i];
    }
    delete [] id_restart_peratom;
    delete [] style_restart_peratom;
    delete [] index_restart_peratom;
  }

  nfix_restart_global = nfix_restart_peratom = 0;

  if(0 == n_ms && have_ms_in_restart)
    error->all(FLERR,"Restart data contains multi-sphere data, which was not restarted. In order to restart it,\n"
                         "you have to place a fix multisphere/* command before the first run command in the input script\n");
}

/* ----------------------------------------------------------------------
   create list of fix indices for fixes which match mask
------------------------------------------------------------------------- */

void Modify::list_init(int mask, int &n, int *&list)
{
  delete [] list;

  n = 0;
  for (int i = 0; i < nfix; i++) if (fmask[i] & mask) n++;
  list = new int[n];

  n = 0;
  for (int i = 0; i < nfix; i++) if (fmask[i] & mask) list[n++] = i;
}

/* ----------------------------------------------------------------------
   create list of fix indices for for pre_exchange fixes
   have contacthistory fixes always come first so it can copy the data
------------------------------------------------------------------------- */

void Modify::list_init_pre_exchange(int mask, int &n, int *&list)
{
  delete [] list;

  n = 0;
  for (int i = 0; i < nfix; i++) if (fmask[i] & mask) n++;
  list = new int[n];

  n = 0;

  for (int i = 0; i < nfix; i++) if (fmask[i] & mask)
  {
    if(0 == strncmp(fix[i]->style,"contacthistory",14))
    //if(0 == strcmp(fix[i]->style,"contacthistory"))
        list[n++] = i;
  }

  for (int i = 0; i < nfix; i++)
  {
      if(0 == strncmp(fix[i]->style,"contacthistory",14))
      //if(0 == strcmp(fix[i]->style,"contacthistory"))
        continue;

      if (fmask[i] & mask) list[n++] = i;
  }
}

/* ----------------------------------------------------------------------
   create list of fix indices for end_of_step fixes
   also create end_of_step_every[]
------------------------------------------------------------------------- */

void Modify::list_init_end_of_step(int mask, int &n, int *&list)
{
  delete [] list;
  delete [] end_of_step_every;

  n = 0;
  for (int i = 0; i < nfix; i++) if (fmask[i] & mask) n++;
  list = new int[n];
  end_of_step_every = new int[n];

  n = 0;
  for (int i = 0; i < nfix; i++)
    if (fmask[i] & mask) {
      list[n] = i;
      end_of_step_every[n++] = fix[i]->nevery;
    }
}

/* ----------------------------------------------------------------------
   create list of fix indices for thermo energy fixes
   only added to list if fix has THERMO_ENERGY mask
   and its thermo_energy flag was set via fix_modify
------------------------------------------------------------------------- */

void Modify::list_init_thermo_energy(int mask, int &n, int *&list)
{
  delete [] list;

  n = 0;
  for (int i = 0; i < nfix; i++)
    if (fmask[i] & mask && fix[i]->thermo_energy) n++;
  list = new int[n];

  n = 0;
  for (int i = 0; i < nfix; i++)
    if (fmask[i] & mask && fix[i]->thermo_energy) list[n++] = i;
}

/* ----------------------------------------------------------------------
   create list of compute indices for computes which store invocation times
------------------------------------------------------------------------- */

void Modify::list_init_compute()
{
  delete [] list_timeflag;

  n_timeflag = 0;
  for (int i = 0; i < ncompute; i++)
    if (compute[i]->timeflag) n_timeflag++;
  list_timeflag = new int[n_timeflag];

  n_timeflag = 0;
  for (int i = 0; i < ncompute; i++)
    if (compute[i]->timeflag) list_timeflag[n_timeflag++] = i;
}

/* ----------------------------------------------------------------------
   return # of bytes of allocated memory from all fixes
------------------------------------------------------------------------- */

bigint Modify::memory_usage()
{
  bigint bytes = 0;
  for (int i = 0; i < nfix; i++)
    bytes += static_cast<bigint> (fix[i]->memory_usage());
  for (int i = 0; i < ncompute; i++)
    bytes += static_cast<bigint> (compute[i]->memory_usage());
  return bytes;
}

/* ======================================================================
   helper functions by Richard Berger (JKU)
========================================================================= */

/* ----------------------------------------------------------------------
   calls a member method on all fixes
------------------------------------------------------------------------- */

void Modify::call_method_on_fixes(FixMethod method) {
  if(timing) {
    for (int i = 0; i < nfix; i++) {
      fix[i]->begin_time_recording();
      (fix[i]->*method)();
      fix[i]->end_time_recording();
    }
  }
  else
  {
    for (int i = 0; i < nfix; i++) {
      (fix[i]->*method)();
    }
  }
}

/* ----------------------------------------------------------------------
   calls a member method on all fixes in the specified list
------------------------------------------------------------------------- */

void Modify::call_method_on_fixes(FixMethod method, int *& ilist, int & inum) {
  if(timing) {
    for (int i = 0; i < inum; i++) {
      const int ifix = ilist[i];
      fix[ifix]->begin_time_recording();
      (fix[ifix]->*method)();
      fix[ifix]->end_time_recording();
    }
  }
  else
  {
    for (int i = 0; i < inum; i++) {
      (fix[ilist[i]]->*method)();
    }
  }
}

/* ----------------------------------------------------------------------
   calls a member method with vflag parameter on all fixes
------------------------------------------------------------------------- */

void Modify::call_method_on_fixes(FixMethodWithVFlag method, int vflag) {
  if(timing) {
    for (int i = 0; i < nfix; i++) {
      fix[i]->begin_time_recording();
      (fix[i]->*method)(vflag);
      fix[i]->end_time_recording();
    }
  }
  else
  {
    for (int i = 0; i < nfix; i++) {
      (fix[i]->*method)(vflag);
    }
  }
}

/* ----------------------------------------------------------------------
   calls a member method with vflag parameter on all fixes in the
   specified list
------------------------------------------------------------------------- */

void Modify::call_method_on_fixes(FixMethodWithVFlag method, int vflag, int *& ilist, int & inum) {
  if(timing) {
    for (int i = 0; i < inum; i++) {
      const int ifix = ilist[i];
      fix[ifix]->begin_time_recording();
      (fix[ifix]->*method)(vflag);
      fix[ifix]->end_time_recording();
    }
  }
  else
  {
    for (int i = 0; i < inum; i++) {
      (fix[ilist[i]]->*method)(vflag);
    }
  }
}

/* ----------------------------------------------------------------------
   calls a respa member method with 2 int parameters on all fixes in the
   specified list
------------------------------------------------------------------------- */

void Modify::call_respa_method_on_fixes(FixMethodRESPA2 method,
    int arg1, int arg2, int *& ilist, int & inum) {
  if(timing) {
    for (int i = 0; i < inum; i++) {
      const int ifix = ilist[i];
      fix[ifix]->begin_time_recording();
      (fix[ifix]->*method)(arg1, arg2);
      fix[ifix]->end_time_recording();
    }
  }
  else
  {
    for (int i = 0; i < inum; i++) {
      (fix[ilist[i]]->*method)(arg1, arg2);
    }
  }
}

/* ----------------------------------------------------------------------
   calls a respa member method with 3 int parameters on all fixes in the
   specified list
------------------------------------------------------------------------- */

void Modify::call_respa_method_on_fixes(FixMethodRESPA3 method, int arg1,
    int arg2, int arg3, int *& ilist, int & inum) {
  if(timing) {
    for (int i = 0; i < inum; i++) {
      const int ifix = ilist[i];
      fix[ifix]->begin_time_recording();
      (fix[ifix]->*method)(arg1, arg2, arg3);
      fix[ifix]->end_time_recording();
    }
  }
  else
  {
    for (int i = 0; i < inum; i++) {
      (fix[ilist[i]]->*method)(arg1, arg2, arg3);
    }
  }
}

/* ----------------------------------------------------------------------
   Updates all computes that requested it at the end of a run
------------------------------------------------------------------------- */

void Modify::update_computes_on_run_end()
{
    for (int i = 0; i < ncompute; i++)
    {
        if (compute[i]->update_on_run_end())
        {
            if (compute[i]->scalar_flag)
            {
                if (!(compute[i]->invoked_flag & INVOKED_SCALAR))
                {
                  compute[i]->compute_scalar();
                  compute[i]->invoked_flag |= INVOKED_SCALAR;
                }
            }
            if (compute[i]->vector_flag)
            {
                if (!(compute[i]->invoked_flag & INVOKED_VECTOR))
                {
                  compute[i]->compute_vector();
                  compute[i]->invoked_flag |= INVOKED_VECTOR;
                }
            }
            if (compute[i]->array_flag)
            {
                if (!(compute[i]->invoked_flag & INVOKED_ARRAY))
                {
                  compute[i]->compute_array();
                  compute[i]->invoked_flag |= INVOKED_ARRAY;
                }
            }
        }
    }
}
