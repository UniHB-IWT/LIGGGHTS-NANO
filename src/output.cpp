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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "output.h"
#include "style_dump.h"
#include "atom.h"
#include "neighbor.h"
#include "input.h"
#include "variable.h"
#include "comm.h"
#include "update.h"
#include "group.h"
#include "domain.h"
#include "thermo.h"
#include "modify.h"
#include "compute.h"
#include "force.h"
#include "dump.h"
#include "write_restart.h"
#include "accelerator_cuda.h"
#include "memory.h"
#include "error.h"
#include "signal_handling.h"

using namespace LAMMPS_NS;

#define DELTA 1

/* ----------------------------------------------------------------------
   initialize all output
------------------------------------------------------------------------- */

Output::Output(LAMMPS *lmp) :
    Pointers(lmp),
    restart_flag(false),
    restart_flag_single(false),
    restart_flag_double(false),
    next_restart(0),
    next_restart_single(0),
    next_restart_double(0),
    restart_every_single(0),
    restart_every_double(0),
    last_restart(0),
    restart_toggle(0),
    var_restart_single(NULL),
    var_restart_double(NULL),
    ivar_restart_single(0),
    ivar_restart_double(0),
    restart1(NULL),
    restart2a(NULL),
    restart2b(NULL),
    restart(NULL)
{
  char **newarg = new char*[4];
  // create a default compute that calculates the temperature of the system
  // NOTE: This compute is deprecated and will be removed in the future
  newarg[0] = (char *) "thermo_temp";
  newarg[1] = (char *) "all";
  newarg[2] = (char *) "temp";
  modify->add_compute(3,newarg,lmp->suffix);
  // create a default compute that calculates the kinetic energy of the system
  newarg[0] = (char *) "thermo_kin_eng";
  newarg[1] = (char *) "all";
  newarg[2] = (char *) "ke";
  modify->add_compute(3,newarg,lmp->suffix);
  delete [] newarg;

  // create default Thermo class
  newarg = new char*[1];
  newarg[0] = (char *) "one";
  thermo = new Thermo(lmp,1,newarg);
  delete [] newarg;

  thermo_every = 0;
  var_thermo = NULL;

  ndump = 0;
  max_dump = 0;
  every_dump = NULL;
  next_dump = NULL;
  last_dump = NULL;
  var_dump = NULL;
  ivar_dump = NULL;
  dump = NULL;
}

/* ----------------------------------------------------------------------
   free all memory
------------------------------------------------------------------------- */

Output::~Output()
{
    if (thermo)
        delete thermo;
    delete [] var_thermo;

    memory->destroy(every_dump);
    memory->destroy(next_dump);
    memory->destroy(last_dump);
    for (int i = 0; i < ndump; i++)
        delete [] var_dump[i];
    memory->sfree(var_dump);
    memory->destroy(ivar_dump);
    for (int i = 0; i < ndump; i++)
        delete dump[i];
    memory->sfree(dump);

    delete [] restart1;
    delete [] restart2a;
    delete [] restart2b;
    delete [] var_restart_single;
    delete [] var_restart_double;
    delete restart;
}

/* ---------------------------------------------------------------------- */

void Output::init()
{
  thermo->init();
  if (var_thermo) {
    ivar_thermo = input->variable->find(var_thermo);
    if (ivar_thermo < 0)
      error->all(FLERR,"Variable name for thermo every does not exist");
    if (!input->variable->equalstyle(ivar_thermo))
      error->all(FLERR,"Variable for thermo every is invalid style");
  }

  for (int i = 0; i < ndump; i++) dump[i]->init();
  for (int i = 0; i < ndump; i++)
    if (every_dump[i] == 0) {
      ivar_dump[i] = input->variable->find(var_dump[i]);
      if (ivar_dump[i] < 0)
        error->all(FLERR,"Variable name for dump every does not exist");
      if (!input->variable->equalstyle(ivar_dump[i]))
        error->all(FLERR,"Variable for dump every is invalid style");
    }

  if (restart_flag_single && restart_every_single == 0) {
    ivar_restart_single = input->variable->find(var_restart_single);
    if (ivar_restart_single < 0)
      error->all(FLERR,"Variable name for restart does not exist");
    if (!input->variable->equalstyle(ivar_restart_single))
      error->all(FLERR,"Variable for restart is invalid style");
  }
  if (restart_flag_double && restart_every_double == 0) {
    ivar_restart_double = input->variable->find(var_restart_double);
    if (ivar_restart_double < 0)
      error->all(FLERR,"Variable name for restart does not exist");
    if (!input->variable->equalstyle(ivar_restart_double))
      error->all(FLERR,"Variable for restart is invalid style");
  }
}

/* ----------------------------------------------------------------------
   perform output for setup of run/min
   do dump first, so memory_usage will include dump allocation
   do thermo last, so will print after memory_usage
   memflag = 0/1 for printing out memory usage
------------------------------------------------------------------------- */

void Output::setup(int memflag)
{
  bigint ntimestep = update->ntimestep;

  // perform dump at start of run only if:
  //   current timestep is multiple of every and last dump not >= this step
  //   this is first run after dump created and firstflag is set
  //   note that variable freq will not write unless triggered by firstflag
  // set next_dump to multiple of every or variable value
  // set next_dump_any to smallest next_dump
  // wrap dumps that invoke computes and variable eval with clear/add
  // if dump not written now, use addstep_compute_all() since don't know
  //   what computes the dump write would invoke
  // if no dumps, set next_dump_any to last+1 so will not influence next

  int writeflag;

  if (ndump && update->restrict_output == 0) {
    for (int idump = 0; idump < ndump; idump++) {
      if (dump[idump]->clearstep || every_dump[idump] == 0)
        modify->clearstep_compute();
      writeflag = 0;
      if (every_dump[idump] && ntimestep % every_dump[idump] == 0 &&
          last_dump[idump] != ntimestep) writeflag = 1;
      if (last_dump[idump] < 0 && dump[idump]->first_flag == 1) writeflag = 1;

      if (writeflag) {
        dump[idump]->write();
        last_dump[idump] = ntimestep;
      }
      if (every_dump[idump])
        next_dump[idump] =
          (ntimestep/every_dump[idump])*every_dump[idump] + every_dump[idump];
      else {
        bigint nextdump = static_cast<bigint>
          (input->variable->compute_equal(ivar_dump[idump]));
        if (nextdump <= ntimestep)
          error->all(FLERR,"Dump every variable returned a bad timestep");
        next_dump[idump] = nextdump;
      }
      if (dump[idump]->clearstep || every_dump[idump] == 0) {
        if (writeflag) modify->addstep_compute(next_dump[idump]);
        else modify->addstep_compute_all(next_dump[idump]);
      }
      if (idump) next_dump_any = MIN(next_dump_any,next_dump[idump]);
      else next_dump_any = next_dump[0];
    }
  } else next_dump_any = update->laststep + 1;

  // do not write restart files at start of run
  // set next_restart values to multiple of every or variable value
  // wrap variable eval with clear/add
  // if no restarts, set next_restart to last+1 so will not influence next

  if (restart_flag && update->restrict_output == 0) {
    if (restart_flag_single) {
      if (restart_every_single)
        next_restart_single =
          (ntimestep/restart_every_single)*restart_every_single +
          restart_every_single;
      else {
        bigint nextrestart = static_cast<bigint>
          (input->variable->compute_equal(ivar_restart_single));
        if (nextrestart <= ntimestep)
          error->all(FLERR,"Restart variable returned a bad timestep");
        next_restart_single = nextrestart;
      }
    } else next_restart_single = update->laststep + 1;
    if (restart_flag_double) {
      if (restart_every_double)
        next_restart_double =
          (ntimestep/restart_every_double)*restart_every_double +
          restart_every_double;
      else {
        bigint nextrestart = static_cast<bigint>
          (input->variable->compute_equal(ivar_restart_double));
        if (nextrestart <= ntimestep)
          error->all(FLERR,"Restart variable returned a bad timestep");
        next_restart_double = nextrestart;
      }
    } else next_restart_double = update->laststep + 1;
    next_restart = MIN(next_restart_single,next_restart_double);
  } else next_restart = update->laststep + 1;

  // print memory usage unless being called between multiple runs

  if (memflag) memory_usage();

  // set next_thermo to multiple of every or variable eval if var defined
  // insure thermo output on last step of run
  // thermo may invoke computes so wrap with clear/add

  modify->clearstep_compute();

  thermo->header();
  thermo->compute(0);
  last_thermo = ntimestep;

  if (var_thermo) {
    next_thermo = static_cast<bigint>
      (input->variable->compute_equal(ivar_thermo));
    if (next_thermo <= ntimestep)
      error->all(FLERR,"Thermo every variable returned a bad timestep");
  } else if (thermo_every) {
    next_thermo = (ntimestep/thermo_every)*thermo_every + thermo_every;
    next_thermo = MIN(next_thermo,update->laststep);
  } else next_thermo = update->laststep;

  modify->addstep_compute(next_thermo);

  // next = next timestep any output will be done

  next = MIN(next_dump_any,next_restart);
  next = MIN(next,next_thermo);
}

/* ----------------------------------------------------------------------
   perform all output for this timestep
   only perform output if next matches current step and last output doesn't
   do dump/restart before thermo so thermo CPU time will include them
------------------------------------------------------------------------- */

void Output::write(bigint ntimestep)
{
  // next_dump does not force output on last step of run
  // wrap dumps that invoke computes or eval of variable with clear/add
  // download data from GPU if necessary

  if (next_dump_any == ntimestep) {
    if (lmp->cuda && !lmp->cuda->oncpu) lmp->cuda->downloadAll();

    for (int idump = 0; idump < ndump; idump++) {
      if (next_dump[idump] == ntimestep) {
        if (dump[idump]->clearstep || every_dump[idump] == 0)
          modify->clearstep_compute();
        if (last_dump[idump] != ntimestep) {
          dump[idump]->write();
          last_dump[idump] = ntimestep;
        }
        if (every_dump[idump]) next_dump[idump] += every_dump[idump];
        else {
          bigint nextdump = static_cast<bigint>
            (input->variable->compute_equal(ivar_dump[idump]));
          if (nextdump <= ntimestep)
            error->all(FLERR,"Dump every variable returned a bad timestep");
          next_dump[idump] = nextdump;
        }
        if (dump[idump]->clearstep || every_dump[idump] == 0)
          modify->addstep_compute(next_dump[idump]);
      }
      if (idump) next_dump_any = MIN(next_dump_any,next_dump[idump]);
      else next_dump_any = next_dump[0];
    }
  }

  // next_restart does not force output on last step of run
  // for toggle = 0, replace "*" with current timestep in restart filename
  // download data from GPU if necessary
  // eval of variable may invoke computes so wrap with clear/add

  if (next_restart == ntimestep) {
    if (lmp->cuda && !lmp->cuda->oncpu) lmp->cuda->downloadAll();

    if (next_restart_single == ntimestep) {
      char *file = new char[strlen(restart1) + 16];
      char *ptr = strchr(restart1,'*');
      *ptr = '\0';
      sprintf(file,"%s" BIGINT_FORMAT "%s",restart1,ntimestep,ptr+1);
      *ptr = '*';
      if (last_restart != ntimestep) restart->write(file);
      delete [] file;
      if (restart_every_single) next_restart_single += restart_every_single;
      else {
        modify->clearstep_compute();
        bigint nextrestart = static_cast<bigint>
          (input->variable->compute_equal(ivar_restart_single));
        if (nextrestart <= ntimestep)
          error->all(FLERR,"Restart variable returned a bad timestep");
        next_restart_single = nextrestart;
        modify->addstep_compute(next_restart_single);
      }
    }
    if (next_restart_double == ntimestep) {
      if (last_restart != ntimestep) {
        if (restart_toggle == 0) {
          restart->write(restart2a);
          restart_toggle = 1;
        } else {
          restart->write(restart2b);
          restart_toggle = 0;
        }
      }
      if (restart_every_double) next_restart_double += restart_every_double;
      else {
        modify->clearstep_compute();
        bigint nextrestart = static_cast<bigint>
          (input->variable->compute_equal(ivar_restart_double));
        if (nextrestart <= ntimestep)
          error->all(FLERR,"Restart variable returned a bad timestep");
        next_restart_double = nextrestart;
        modify->addstep_compute(next_restart_double);
      }
    }
    last_restart = ntimestep;
    next_restart = MIN(next_restart_single,next_restart_double);

    if (SignalHandler::request_write_restart) {
        char *file = new char[24 + 16 + 5];
        sprintf(file,"restart_forced_liggghts_" BIGINT_FORMAT ".data",ntimestep);
        bool has_restart = restart != NULL;
        if (!has_restart)
            restart = new WriteRestart(lmp);
        restart->write(file);
        if (!has_restart)
        {
            delete restart;
            restart = NULL;
        }
        delete [] file;
        SignalHandler::request_write_restart = false;
        error->warning(FLERR, "Forced restart written");
    }
  }

  // insure next_thermo forces output on last step of run
  // thermo may invoke computes so wrap with clear/add

  if (next_thermo == ntimestep) {
    modify->clearstep_compute();
    // check all computes and those with update_on_run_end activated will be updated
    if (ntimestep == update->laststep)
        modify->update_computes_on_run_end();
    if (last_thermo != ntimestep) thermo->compute(1);
    last_thermo = ntimestep;
    if (var_thermo) {
      next_thermo = static_cast<bigint>
        (input->variable->compute_equal(ivar_thermo));
      if (next_thermo <= ntimestep)
        error->all(FLERR,"Thermo every variable returned a bad timestep");
    } else if (thermo_every) next_thermo += thermo_every;
    else next_thermo = update->laststep;
    next_thermo = MIN(next_thermo,update->laststep);
    modify->addstep_compute(next_thermo);
  }

  // next = next timestep any output will be done

  next = MIN(next_dump_any,next_restart);
  next = MIN(next,next_thermo);
}

/* ----------------------------------------------------------------------
   force a snapshot to be written for all dumps
   called from PRD and TAD
------------------------------------------------------------------------- */

void Output::write_dump(bigint ntimestep)
{
  for (int idump = 0; idump < ndump; idump++) {
    dump[idump]->write();
    last_dump[idump] = ntimestep;
  }
}

/* ----------------------------------------------------------------------
   force restart file(s) to be written
   called from PRD and TAD
------------------------------------------------------------------------- */

void Output::write_restart(bigint ntimestep)
{
  if (restart_flag_single) {
    char *file = new char[strlen(restart1) + 16];
    char *ptr = strchr(restart1,'*');
    *ptr = '\0';
    sprintf(file,"%s" BIGINT_FORMAT "%s",restart1,ntimestep,ptr+1);
    *ptr = '*';
    restart->write(file);
    delete [] file;
  }

  if (restart_flag_double) {
    if (restart_toggle == 0) {
      restart->write(restart2a);
      restart_toggle = 1;
    } else {
      restart->write(restart2b);
      restart_toggle = 0;
    }
  }

  last_restart = ntimestep;
}

/* ----------------------------------------------------------------------
   timestep is being changed, called by update->reset_timestep()
   reset next timestep values for dumps, restart, thermo output
   reset to smallest value >= new timestep
   if next timestep set by variable evaluation,
     eval for ntimestep-1, so current ntimestep can be returned if needed
     no guarantee that variable can be evaluated for ntimestep-1
       if it depends on computes, but live with that rare case for now
------------------------------------------------------------------------- */

void Output::reset_timestep(bigint ntimestep)
{
  next_dump_any = MAXBIGINT;
  for (int idump = 0; idump < ndump; idump++) {
    if (every_dump[idump]) {
      next_dump[idump] = (ntimestep/every_dump[idump])*every_dump[idump];
      if (next_dump[idump] < ntimestep) next_dump[idump] += every_dump[idump];
    } else {
      modify->clearstep_compute();
      update->ntimestep--;
      bigint nextdump = static_cast<bigint>
        (input->variable->compute_equal(ivar_dump[idump]));
      if (nextdump < ntimestep)
        error->all(FLERR,"Dump every variable returned a bad timestep");
      update->ntimestep++;
      next_dump[idump] = nextdump;
      modify->addstep_compute(next_dump[idump]);
    }
    next_dump_any = MIN(next_dump_any,next_dump[idump]);
  }

  if (restart_flag_single) {
    if (restart_every_single) {
      next_restart_single =
        (ntimestep/restart_every_single)*restart_every_single;
      if (next_restart_single < ntimestep)
        next_restart_single += restart_every_single;
    } else {
      modify->clearstep_compute();
      update->ntimestep--;
      bigint nextrestart = static_cast<bigint>
        (input->variable->compute_equal(ivar_restart_single));
      if (nextrestart < ntimestep)
        error->all(FLERR,"Restart variable returned a bad timestep");
      update->ntimestep++;
      next_restart_single = nextrestart;
      modify->addstep_compute(next_restart_single);
    }
  } else next_restart_single = update->laststep + 1;

  if (restart_flag_double) {
    if (restart_every_double) {
      next_restart_double =
        (ntimestep/restart_every_double)*restart_every_double;
      if (next_restart_double < ntimestep)
        next_restart_double += restart_every_double;
    } else {
      modify->clearstep_compute();
      update->ntimestep--;
      bigint nextrestart = static_cast<bigint>
        (input->variable->compute_equal(ivar_restart_double));
      if (nextrestart < ntimestep)
        error->all(FLERR,"Restart variable returned a bad timestep");
      update->ntimestep++;
      next_restart_double = nextrestart;
      modify->addstep_compute(next_restart_double);
    }
  } else next_restart_double = update->laststep + 1;

  next_restart = MIN(next_restart_single,next_restart_double);

  if (var_thermo) {
    modify->clearstep_compute();
    update->ntimestep--;
    next_thermo = static_cast<bigint>
      (input->variable->compute_equal(ivar_thermo));
    if (next_thermo < ntimestep)
      error->all(FLERR,"Thermo every variable returned a bad timestep");
    update->ntimestep++;
    next_thermo = MIN(next_thermo,update->laststep);
    modify->addstep_compute(next_thermo);
  } else if (thermo_every) {
    next_thermo = (ntimestep/thermo_every)*thermo_every;
    if (next_thermo < ntimestep) next_thermo += thermo_every;
    next_thermo = MIN(next_thermo,update->laststep);
  } else next_thermo = update->laststep;

  next = MIN(next_dump_any,next_restart);
  next = MIN(next,next_thermo);
}

/* ----------------------------------------------------------------------
   add a Dump to list of Dumps
------------------------------------------------------------------------- */

void Output::add_dump(int narg, char **arg)
{
  if (narg < 5) error->all(FLERR,"Illegal dump command");

  // error checks

  for (int idump = 0; idump < ndump; idump++)
    if (strcmp(arg[0],dump[idump]->id) == 0)
      error->all(FLERR,"Reuse of dump ID");
  int igroup = group->find(arg[1]);
  if (igroup == -1) error->all(FLERR,"Could not find dump group ID");
  if (force->inumeric(FLERR,arg[3]) <= 0)
    error->all(FLERR,"Invalid dump frequency");

  // extend Dump list if necessary

  if (ndump == max_dump) {
    max_dump += DELTA;
    dump = (Dump **)
      memory->srealloc(dump,max_dump*sizeof(Dump *),"output:dump");
    memory->grow(every_dump,max_dump,"output:every_dump");
    memory->grow(next_dump,max_dump,"output:next_dump");
    memory->grow(last_dump,max_dump,"output:last_dump");
    var_dump = (char **)
      memory->srealloc(var_dump,max_dump*sizeof(char *),"output:var_dump");
    memory->grow(ivar_dump,max_dump,"output:ivar_dump");
  }

  // create the Dump

  if (0) return;         // dummy line to enable else-if macro expansion

#define DUMP_CLASS
#define DumpStyle(key,Class) \
  else if (strcmp(arg[2],#key) == 0) dump[ndump] = new Class(lmp,narg,arg);
#include "style_dump.h"
#undef DUMP_CLASS

  else error->all(FLERR,"Invalid dump style");

  every_dump[ndump] = force->inumeric(FLERR,arg[3]);
  if (every_dump[ndump] <= 0) error->all(FLERR,"Illegal dump command");
  last_dump[ndump] = -1;
  var_dump[ndump] = NULL;
  ndump++;
}

/* ----------------------------------------------------------------------
   modify parameters of a Dump
------------------------------------------------------------------------- */

void Output::modify_dump(int narg, char **arg)
{
  if (narg < 1) error->all(FLERR,"Illegal dump_modify command");

  // find which dump it is

  int idump;
  for (idump = 0; idump < ndump; idump++)
    if (strcmp(arg[0],dump[idump]->id) == 0) break;
  if (idump == ndump) error->all(FLERR,"Cound not find dump_modify ID");

  dump[idump]->modify_params(narg-1,&arg[1]);
}

/* ----------------------------------------------------------------------
   delete a Dump from list of Dumps
------------------------------------------------------------------------- */

void Output::delete_dump(char *id)
{
  // find which dump it is and delete it

  int idump;
  for (idump = 0; idump < ndump; idump++)
    if (strcmp(id,dump[idump]->id) == 0) break;
  if (idump == ndump) error->all(FLERR,"Could not find undump ID");

  delete dump[idump];
  delete [] var_dump[idump];

  // move other dumps down in list one slot

  for (int i = idump+1; i < ndump; i++) {
    dump[i-1] = dump[i];
    every_dump[i-1] = every_dump[i];
    next_dump[i-1] = next_dump[i];
    last_dump[i-1] = last_dump[i];
    var_dump[i-1] = var_dump[i];
    ivar_dump[i-1] = ivar_dump[i];
  }
  ndump--;
}

/* ----------------------------------------------------------------------
   set thermo output frequency from input script
------------------------------------------------------------------------- */

void Output::set_thermo(int narg, char **arg)
{
  if (narg != 1) error->all(FLERR,"Illegal thermo command");

  if (strstr(arg[0],"v_") == arg[0]) {
    delete [] var_thermo;
    int n = strlen(&arg[0][2]) + 1;
    var_thermo = new char[n];
    strcpy(var_thermo,&arg[0][2]);
  } else {
    thermo_every = force->inumeric(FLERR,arg[0]);
    if (thermo_every < 0) error->all(FLERR,"Illegal thermo command");
  }
}

/* ----------------------------------------------------------------------
   new Thermo style
------------------------------------------------------------------------- */

void Output::create_thermo(int narg, char **arg)
{
  if (narg < 1) error->all(FLERR,"Illegal thermo_style command");

  // don't allow this so that dipole style can safely allocate inertia vector

  if (domain->box_exist == 0)
    error->all(FLERR,"Thermo_style command before simulation box is defined");

  // warn if previous thermo had been modified via thermo_modify command

  if (thermo->modified && comm->me == 0 && !lmp->wb)
    error->warning(FLERR,"New thermo_style command, "
                   "previous thermo_modify settings will be lost");

  // set thermo = NULL in case new Thermo throws an error

  delete thermo;
  thermo = NULL;
  thermo = new Thermo(lmp,narg,arg);
}

/* ----------------------------------------------------------------------
   setup restart capability for single or double output files
   if only one filename and it contains no "*", then append ".*"
------------------------------------------------------------------------- */

void Output::create_restart(int narg, char **arg)
{
  if (narg < 1) error->all(FLERR,"Illegal restart command");

  int every = 0;
  int varflag = 0;

  if (strstr(arg[0],"v_") == arg[0]) varflag = 1;
  else every = force->inumeric(FLERR,arg[0]);

  if (!varflag && every == 0) {
    if (narg != 1) error->all(FLERR,"Illegal restart command");

    restart_flag = restart_flag_single = restart_flag_double = false;
    last_restart = -1;

    delete restart;
    restart = NULL;
    delete [] restart1;
    delete [] restart2a;
    delete [] restart2b;
    restart1 = restart2a = restart2b = NULL;
    delete [] var_restart_single;
    delete [] var_restart_double;
    var_restart_single = var_restart_double = NULL;

    return;
  }

  if (narg != 2 && narg != 3) error->all(FLERR,"Illegal restart command");

  if (narg == 2) {
    restart_flag = restart_flag_single = true;

    if (varflag) {
      delete [] var_restart_single;
      int n = strlen(&arg[0][2]) + 1;
      var_restart_single = new char[n];
      strcpy(var_restart_single,&arg[0][2]);
      restart_every_single = 0;
    } else restart_every_single = every;

    int n = strlen(arg[1]) + 3;
    restart1 = new char[n];
    strcpy(restart1,arg[1]);
    if (strchr(restart1,'*') == NULL) strcat(restart1,".*");
  }

  if (narg == 3) {
    restart_flag = restart_flag_double = true;

    if (varflag) {
      delete [] var_restart_double;
      int n = strlen(&arg[0][2]) + 1;
      var_restart_double = new char[n];
      strcpy(var_restart_double,&arg[0][2]);
      restart_every_double = 0;
    } else restart_every_double = every;

    restart_toggle = 0;
    int n = strlen(arg[1]) + 3;
    restart2a = new char[n];
    strcpy(restart2a,arg[1]);
    n = strlen(arg[2]) + 1;
    restart2b = new char[n];
    strcpy(restart2b,arg[2]);
  }

  if (restart == NULL) restart = new WriteRestart(lmp);
}

/* ----------------------------------------------------------------------
   sum and print memory usage
   result is only memory on proc 0, not averaged across procs
------------------------------------------------------------------------- */

void Output::memory_usage()
{
  bigint bytes = 0;
  bytes += atom->memory_usage();
  bytes += neighbor->memory_usage();
  bytes += comm->memory_usage();
  bytes += update->memory_usage();
  bytes += force->memory_usage();
  bytes += modify->memory_usage();
  for (int i = 0; i < ndump; i++) bytes += dump[i]->memory_usage();

  double mbytes = bytes/1024.0/1024.0;

  if (comm->me == 0) {
    if (screen)
      fprintf(screen,"Memory usage per processor = %g Mbytes\n",mbytes);
    if (logfile)
      fprintf(logfile,"Memory usage per processor = %g Mbytes\n",mbytes);
  }
}

/* ----------------------------------------------------------------------
   identifies when the next restart will be written
   also handles signals to force writing of a restart
   this function is called by Neighbor::decide
------------------------------------------------------------------------- */

bool Output::restart_requested(const bigint ntimestep)
{
    if (SignalHandler::request_write_restart || SignalHandler::request_quit)
    {
        // we have something to write now
        next = ntimestep;
        // if quit is request write thermo
        if (SignalHandler::request_quit)
            next_thermo = ntimestep;
        // if restart writing is request do it
        if (SignalHandler::request_write_restart)
            next_restart = ntimestep;
    }
    return next_restart == ntimestep;
}

/* ----------------------------------------------------------------------
   request a restart for a certain timestep
------------------------------------------------------------------------- */

void Output::request_restart(const bigint ntimestep)
{
    if (restart_flag)
    {
        next_restart = ntimestep;
        if (restart_every_single)
            next_restart_single = ntimestep;
        if (restart_every_double)
            next_restart_double = ntimestep;
    }
}
