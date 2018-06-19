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

/* ----------------------------------------------------------------------
   Contributing author: Pieter in 't Veld (SNL)
------------------------------------------------------------------------- */

#include <stdlib.h>
#include <string.h>
#include "unistd.h"
#include "fix_ave_spatial.h"
#include "atom.h"
#include "update.h"
#include "force.h"
#include "domain.h"
#include "region.h"
#include "lattice.h"
#include "modify.h"
#include "compute.h"
#include "input.h"
#include "variable.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace FixConst;

enum{LOWER,CENTER,UPPER,COORD};
enum{V,F,DENSITY_NUMBER,DENSITY_MASS,COMPUTE,FIX,VARIABLE};
enum{SAMPLE,ALL};
enum{BOX,LATTICE,REDUCED};
enum{ONE,RUNNING,WINDOW};

#define BIG 1000000000

/* ---------------------------------------------------------------------- */

FixAveSpatial::FixAveSpatial(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg < 6) error->all(FLERR,"Illegal fix ave/spatial command");

  MPI_Comm_rank(world,&me);

  nevery = force->inumeric(FLERR,arg[3]);
  nrepeat = force->inumeric(FLERR,arg[4]);
  nfreq = force->inumeric(FLERR,arg[5]);

  global_freq = nfreq;
  no_change_box = 1;

  write_ts_ = true;

  ndim = 0;
  int iarg = 6;
  while (iarg < narg && ndim < 3) {
    if (iarg+3 > narg) break;
    if (strcmp(arg[iarg],"x") == 0) dim[ndim] = 0;
    else if (strcmp(arg[iarg],"y") == 0) dim[ndim] = 1;
    else if (strcmp(arg[iarg],"z") == 0) dim[ndim] = 2;
    else break;

    if (dim[ndim] == 2 && domain->dimension == 2)
      error->all(FLERR,"Cannot use fix ave/spatial z for 2 dimensional model");

    if (strcmp(arg[iarg+1],"lower") == 0) originflag[ndim] = LOWER;
    else if (strcmp(arg[iarg+1],"center") == 0) originflag[ndim] = CENTER;
    else if (strcmp(arg[iarg+1],"upper") == 0) originflag[ndim] = UPPER;
    else originflag[ndim] = COORD;
    if (originflag[ndim] == COORD)
      origin[ndim] = force->numeric(FLERR,arg[iarg+1]);

    delta[ndim] = force->numeric(FLERR,arg[iarg+2]);
    ndim++;
    iarg += 3;
  }

  if (!ndim) error->all(FLERR,"Illegal fix ave/spatial command");
  if (ndim == 2 && dim[0] == dim[1])
    error->all(FLERR,"Same dimension twice in fix ave/spatial");
  if (ndim == 3 && (dim[0] == dim[1] || dim[1] == dim[2] || dim[0] == dim[2]))
    error->all(FLERR,"Same dimension twice in fix ave/spatial");

  // parse values until one isn't recognized

  which = new int[narg-9];
  argindex = new int[narg-9];
  ids = new char*[narg-9];
  value2index = new int[narg-9];
  nvalues = 0;

  while (iarg < narg) {
    ids[nvalues] = NULL;

    if (strcmp(arg[iarg],"vx") == 0) {
      which[nvalues] = V;
      argindex[nvalues++] = 0;
    } else if (strcmp(arg[iarg],"vy") == 0) {
      which[nvalues] = V;
      argindex[nvalues++] = 1;
    } else if (strcmp(arg[iarg],"vz") == 0) {
      which[nvalues] = V;
      argindex[nvalues++] = 2;

    } else if (strcmp(arg[iarg],"fx") == 0) {
      which[nvalues] = F;
      argindex[nvalues++] = 0;
    } else if (strcmp(arg[iarg],"fy") == 0) {
      which[nvalues] = F;
      argindex[nvalues++] = 1;
    } else if (strcmp(arg[iarg],"fz") == 0) {
      which[nvalues] = F;
      argindex[nvalues++] = 2;

    } else if (strcmp(arg[iarg],"density/number") == 0) {
      which[nvalues] = DENSITY_NUMBER;
      argindex[nvalues++] = 0;
    } else if (strcmp(arg[iarg],"density/mass") == 0) {
      which[nvalues] = DENSITY_MASS;
      argindex[nvalues++] = 0;

    } else if (strncmp(arg[iarg],"c_",2) == 0 ||
               strncmp(arg[iarg],"f_",2) == 0 ||
               strncmp(arg[iarg],"v_",2) == 0) {
      if (arg[iarg][0] == 'c') which[nvalues] = COMPUTE;
      else if (arg[iarg][0] == 'f') which[nvalues] = FIX;
      else if (arg[iarg][0] == 'v') which[nvalues] = VARIABLE;

      int n = strlen(arg[iarg]);
      char *suffix = new char[n];
      strcpy(suffix,&arg[iarg][2]);

      char *ptr = strchr(suffix,'[');
      if (ptr) {
        if (suffix[strlen(suffix)-1] != ']')
          error->all(FLERR,"Illegal fix ave/spatial command");
        argindex[nvalues] = atoi(ptr+1);
        *ptr = '\0';
      } else argindex[nvalues] = 0;

      n = strlen(suffix) + 1;
      ids[nvalues] = new char[n];
      strcpy(ids[nvalues],suffix);
      nvalues++;
      delete [] suffix;

    } else break;

    iarg++;
  }

  // optional args

  normflag = ALL;
  scaleflag = BOX; 
  regionflag = 0;
  idregion = NULL;
  fp = NULL;
  
  fp2 = NULL;
  lowerLimit = 0;
  upperLimit = 0;
  calcStd = 0;
  ave = ONE;
  nwindow = 0;
  overwrite = 0;
  char *title1 = NULL;
  char *title2 = NULL;
  char *title3 = NULL;

  while (iarg < narg) {
    if (strcmp(arg[iarg],"norm") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix ave/spatial command");
      if (strcmp(arg[iarg+1],"all") == 0) normflag = ALL;
      else if (strcmp(arg[iarg+1],"sample") == 0) normflag = SAMPLE;
      else error->all(FLERR,"Illegal fix ave/spatial command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"units") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix ave/spatial command");
      if (strcmp(arg[iarg+1],"box") == 0) scaleflag = BOX;
      else if (strcmp(arg[iarg+1],"lattice") == 0) scaleflag = LATTICE;
      else if (strcmp(arg[iarg+1],"reduced") == 0) scaleflag = REDUCED;
      else error->all(FLERR,"Illegal fix ave/spatial command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"region") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix ave/spatial command");
      iregion = domain->find_region(arg[iarg+1]);
      if (iregion == -1)
        error->all(FLERR,"Region ID for fix ave/spatial does not exist");
      int n = strlen(arg[iarg+1]) + 1;
      idregion = new char[n];
      strcpy(idregion,arg[iarg+1]);
      regionflag = 1;
      iarg += 2;
    } else if (strcmp(arg[iarg],"file") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix ave/spatial command");
      if (me == 0) {
        fp = fopen(arg[iarg+1],"w");
        if (fp == NULL) {
          char str[512];
          sprintf(str,"Cannot open fix ave/spatial file %s",arg[iarg+1]);
          error->one(FLERR,str);
        }
      }
      iarg += 2;
    } else if (strcmp(arg[iarg],"write_ts") == 0) { 
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix ave/spatial command");
      if(strcmp(arg[iarg+1],"yes") == 0)
        write_ts_ = true;
      else if(strcmp(arg[iarg+1],"no") == 0)
        write_ts_ = false;
      else error->all(FLERR,"Illegal fix ave/spatial command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"ave") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix ave/spatial command");
      if (strcmp(arg[iarg+1],"one") == 0) ave = ONE;
      else if (strcmp(arg[iarg+1],"running") == 0) ave = RUNNING;
      else if (strcmp(arg[iarg+1],"window") == 0) ave = WINDOW;
      else error->all(FLERR,"Illegal fix ave/spatial command");
      if (ave == WINDOW) {
        if (iarg+3 > narg) error->all(FLERR,"Illegal fix ave/spatial command");
        nwindow = force->inumeric(FLERR,arg[iarg+2]);
        if (nwindow <= 0) error->all(FLERR,"Illegal fix ave/spatial command");
      }
      iarg += 2;
      if (ave == WINDOW) iarg++;
    } else if (strcmp(arg[iarg],"overwrite") == 0) {
      overwrite = 1;
      iarg += 1;
    } else if (strcmp(arg[iarg],"title1") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix ave/spatial command");
      delete [] title1;
      int n = strlen(arg[iarg+1]) + 1;
      title1 = new char[n];
      strcpy(title1,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"title2") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix ave/spatial command");
      delete [] title2;
      int n = strlen(arg[iarg+1]) + 1;
      title2 = new char[n];
      strcpy(title2,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"title3") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix ave/spatial command");
      delete [] title3;
      int n = strlen(arg[iarg+1]) + 1;
      title3 = new char[n];
      strcpy(title3,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"std") == 0) { 
      if (iarg+4 > narg) error->all(FLERR,"Illegal fix ave/spatial command");
      calcStd = 1;
      lowerLimit = atoi(arg[iarg+1]);
      upperLimit = atoi(arg[iarg+2]);
      if (me == 0) {
        fp2 = fopen(arg[iarg+3],"w");
        if (fp2 == NULL) {
          char str[512];
          sprintf(str,"Cannot open fix ave/spatial file %s",arg[iarg+3]);
          error->one(FLERR,str);
        }
      }
      iarg += 4;
    } else error->all(FLERR,"Illegal fix ave/spatial command");
  }

  // setup and error check

  if (nevery <= 0 || nrepeat <= 0 || nfreq <= 0)
    error->all(FLERR,"Illegal fix ave/spatial command");
  if (nfreq % nevery || (nrepeat-1)*nevery >= nfreq)
    error->all(FLERR,"Illegal fix ave/spatial command");
  if (delta[0] <= 0.0) error->all(FLERR,"Illegal fix ave/spatial command");
  if (ndim >= 2 && delta[1] <= 0.0)
    error->all(FLERR,"Illegal fix ave/spatial command");
  if (ndim == 3 && delta[2] <= 0.0)
    error->all(FLERR,"Illegal fix ave/spatial command");
  if (ave != RUNNING && overwrite)
    error->all(FLERR,"Illegal fix ave/spatial command");

  for (int i = 0; i < nvalues; i++) {
    if (which[i] == COMPUTE) {
      int icompute = modify->find_compute(ids[i]);
      if (icompute < 0)
        error->all(FLERR,"Compute ID for fix ave/spatial does not exist");
      if (modify->compute[icompute]->peratom_flag == 0)
        error->all(FLERR,"Fix ave/spatial compute does not "
                   "calculate per-atom values");
      if (argindex[i] == 0 &&
          modify->compute[icompute]->size_peratom_cols != 0)
        error->all(FLERR,"Fix ave/spatial compute does not "
                   "calculate a per-atom vector");
      if (argindex[i] && modify->compute[icompute]->size_peratom_cols == 0)
        error->all(FLERR,"Fix ave/spatial compute does not "
                   "calculate a per-atom array");
      if (argindex[i] &&
          argindex[i] > modify->compute[icompute]->size_peratom_cols)
        error->all(FLERR,
                   "Fix ave/spatial compute vector is accessed out-of-range");

    } else if (which[i] == FIX) {
      int ifix = modify->find_fix(ids[i]);
      if (ifix < 0)
        error->all(FLERR,"Fix ID for fix ave/spatial does not exist");
      if (modify->fix[ifix]->peratom_flag == 0)
        error->all(FLERR,
                   "Fix ave/spatial fix does not calculate per-atom values");
      if (argindex[i] && modify->fix[ifix]->size_peratom_cols != 0)
        error->all(FLERR,
                   "Fix ave/spatial fix does not calculate a per-atom vector");
      if (argindex[i] && modify->fix[ifix]->size_peratom_cols == 0)
        error->all(FLERR,
                   "Fix ave/spatial fix does not calculate a per-atom array");
      if (argindex[i] && argindex[i] > modify->fix[ifix]->size_peratom_cols)
        error->all(FLERR,"Fix ave/spatial fix vector is accessed out-of-range");
    } else if (which[i] == VARIABLE) {
      int ivariable = input->variable->find(ids[i]);
      if (ivariable < 0)
        error->all(FLERR,"Variable name for fix ave/spatial does not exist");
      if (input->variable->atomstyle(ivariable) == 0)
        error->all(FLERR,"Fix ave/spatial variable is not atom-style variable");
    }
  }

  // print file comment lines

  if (fp && me == 0) {
    if (title1 && strlen (title1) > 1) fprintf(fp,"%s\n",title1);
    else if(!title1) fprintf(fp,"# Spatial-averaged data for fix %s and group %s\n",
                 id,arg[1]);
    if (title2 && strlen (title2) > 1) fprintf(fp,"%s\n",title2);
    else if(!title2) fprintf(fp,"# Timestep Number-of-bins\n");
    if (title3 && strlen (title3) > 1) fprintf(fp,"%s\n",title3);
    else if(!title3) {
      if (ndim == 1) fprintf(fp,"# Bin Coord Ncount");
      else if (ndim == 2) fprintf(fp,"# Bin Coord1 Coord2 Ncount");
      else if (ndim == 3) fprintf(fp,"# Bin Coord1 Coord2 Coord3 Ncount");
      for (int i = 0; i < nvalues; i++) fprintf(fp," %s",arg[6+3*ndim+i]);
      fprintf(fp,"\n");
    }
    filepos = ftell(fp);
  }

  if (fp2 && me == 0) { 
    fprintf(fp2,"# Spatial-averaged data for fix %s and group %s\n",
                 id,arg[1]);
    fprintf(fp2,"# Mean and standard deviation for all bins including between N1=%i and N2=%i particles\n",lowerLimit,upperLimit);
    fprintf(fp2,"# Timestep  Natoms  Nbins  MaxAtomsPerBin  NbinsEmpty  Nbins<N1  Nbins>N2  Nsamples  NatomsPerSample  ");

    for (int i = 0; i < nvalues; i++) {
      fprintf(fp2,"{%s: trueMean samplesMean Std}   ",arg[6+3*ndim+i]);
    }
    fprintf(fp2,"\n");
  }

  fflush(fp2);

  delete [] title1;
  delete [] title2;
  delete [] title3;

  // this fix produces a global array

  array_flag = 1;
  size_array_rows = BIG;
  size_array_cols = 1 + ndim + nvalues;
  extarray = 0;

  // setup scaling

  int triclinic = domain->triclinic;
  if (triclinic == 1 && scaleflag != REDUCED)
    error->all(FLERR,
               "Fix ave/spatial for triclinic boxes requires units reduced");

  if (scaleflag == LATTICE) {
    xscale = domain->lattice->xlattice;
    yscale = domain->lattice->ylattice;
    zscale = domain->lattice->zlattice;
  }
  else xscale = yscale = zscale = 1.0;

  // apply scaling factors

  double scale = 1.0;
  for (int idim = 0; idim < ndim; idim++) {
    if (dim[idim] == 0) scale = xscale;
    else if (dim[idim] == 1) scale = yscale;
    else if (dim[idim] == 2) scale = zscale;
    delta[idim] *= scale;
    if (originflag[idim] == COORD) origin[idim] *= scale;
    invdelta[idim] = 1.0/delta[idim];
  }

  // initializations

  irepeat = 0;
  iwindow = window_limit = 0;
  norm = 0;

  maxvar = 0;
  varatom = NULL;

  maxatom = 0;
  bin = NULL;

  nbins = maxbin = 0;
  count_one = count_many = count_sum = count_total = NULL;
  coord = NULL;
  count_list = NULL;
  values_one = values_many = values_sum = values_total = NULL;
  values_list = NULL;

  // nvalid = next step on which end_of_step does something
  // add nvalid to all computes that store invocation times
  // since don't know a priori which are invoked by this fix
  // once in end_of_step() can set timestep for ones actually invoked

  nvalid = nextvalid();
  modify->addstep_compute_all(nvalid);
}

/* ---------------------------------------------------------------------- */

FixAveSpatial::~FixAveSpatial()
{
  delete [] which;
  delete [] argindex;
  for (int i = 0; i < nvalues; i++) delete [] ids[i];
  delete [] ids;
  delete [] value2index;
  delete [] idregion;

  if (fp && me == 0) fclose(fp);
  if (fp2 && me == 0) fclose(fp2); 

  memory->destroy(varatom);
  memory->destroy(bin);

  memory->destroy(count_one);
  memory->destroy(count_many);
  memory->destroy(count_sum);
  memory->destroy(count_total);
  memory->destroy(coord);
  memory->destroy(count_list);
  memory->destroy(values_one);
  memory->destroy(values_many);
  memory->destroy(values_sum);
  memory->destroy(values_total);
  memory->destroy(values_list);
}

/* ---------------------------------------------------------------------- */

int FixAveSpatial::setmask()
{
  int mask = 0;
  mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixAveSpatial::init()
{
  // set and check validity of region

  if (regionflag) {
    iregion = domain->find_region(idregion);
    if (iregion == -1)
      error->all(FLERR,"Region ID for fix ave/spatial does not exist");
    region = domain->regions[iregion];
  }

  // # of bins cannot vary for ave = RUNNING or WINDOW

  if (ave == RUNNING || ave == WINDOW) {
    if (scaleflag != REDUCED && domain->box_change_size)
      error->all(FLERR,
                 "Fix ave/spatial settings invalid with changing box size");
  }

  // set indices and check validity of all computes,fixes,variables
  // check that fix frequency is acceptable

  for (int m = 0; m < nvalues; m++) {
    if (which[m] == COMPUTE) {
      int icompute = modify->find_compute(ids[m]);
      if (icompute < 0)
        error->all(FLERR,"Compute ID for fix ave/spatial does not exist");
      value2index[m] = icompute;

    } else if (which[m] == FIX) {
      int ifix = modify->find_fix(ids[m]);
      if (ifix < 0)
        error->all(FLERR,"Fix ID for fix ave/spatial does not exist");
      value2index[m] = ifix;

      if (nevery % modify->fix[ifix]->peratom_freq)
        error->all(FLERR,
                   "Fix for fix ave/spatial not computed at compatible time");

    } else if (which[m] == VARIABLE) {
      int ivariable = input->variable->find(ids[m]);
      if (ivariable < 0)
        error->all(FLERR,"Variable name for fix ave/spatial does not exist");
      value2index[m] = ivariable;

    } else value2index[m] = -1;
  }

  // need to reset nvalid if nvalid < ntimestep b/c minimize was performed

  if (nvalid < update->ntimestep) {
    irepeat = 0;
    nvalid = nextvalid();
    modify->addstep_compute_all(nvalid);
  }
}

/* ----------------------------------------------------------------------
   setup initial bins
   only does averaging if nvalid = current timestep
------------------------------------------------------------------------- */

void FixAveSpatial::setup(int vflag)
{
  setup_bins();
  end_of_step();
}

/* ---------------------------------------------------------------------- */

void FixAveSpatial::end_of_step()
{
  int i,j,m,n;

  // skip if not step which requires doing something

  bigint ntimestep = update->ntimestep;
  if (ntimestep != nvalid) return;

  // zero out arrays that accumulate over many samples
  // if box changes, first re-setup bins

  if (irepeat == 0) {
    if (domain->box_change) setup_bins();
    for (m = 0; m < nbins; m++) {
      count_many[m] = count_sum[m] = 0.0;
      for (i = 0; i < nvalues; i++) values_many[m][i] = 0.0;
    }
  }

  // zero out arrays for one sample

  for (m = 0; m < nbins; m++) {
    count_one[m] = 0.0;
    for (i = 0; i < nvalues; i++) values_one[m][i] = 0.0;
  }

  // assign each atom to a bin

  double **x = atom->x;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  if (nlocal > maxatom) {
    maxatom = atom->nmax;
    memory->destroy(bin);
    memory->create(bin,maxatom,"ave/spatial:bin");
  }

  if (ndim == 1) atom2bin1d();
  else if (ndim == 2) atom2bin2d();
  else atom2bin3d();

  // perform the computation for one sample
  // accumulate results of attributes,computes,fixes,variables to local copy
  // sum within each bin, only include atoms in fix group
  // compute/fix/variable may invoke computes so wrap with clear/add

  modify->clearstep_compute();

  for (m = 0; m < nvalues; m++) {
    n = value2index[m];
    j = argindex[m];

    // V,F adds velocities,forces to values

    if (which[m] == V || which[m] == F) {
      double **attribute;
      if (which[m] == V) attribute = atom->v;
      else attribute = atom->f;

      if (regionflag == 0) {
        for (i = 0; i < nlocal; i++)
          if (mask[i] & groupbit)
            values_one[bin[i]][m] += attribute[i][j];
      } else {
        for (i = 0; i < nlocal; i++)
          if (mask[i] & groupbit && region->match(x[i][0],x[i][1],x[i][2]))
            values_one[bin[i]][m] += attribute[i][j];
      }

    // DENSITY_NUMBER adds 1 to values

    } else if (which[m] == DENSITY_NUMBER) {

      if (regionflag == 0) {
        for (i = 0; i < nlocal; i++)
          if (mask[i] & groupbit)
            values_one[bin[i]][m] += 1.0;
      } else {
        for (i = 0; i < nlocal; i++)
          if (mask[i] & groupbit && region->match(x[i][0],x[i][1],x[i][2]))
            values_one[bin[i]][m] += 1.0;
      }

    // DENSITY_MASS adds mass to values

    } else if (which[m] == DENSITY_MASS) {
      int *type = atom->type;
      double *mass = atom->mass;
      double *rmass = atom->rmass;

      if (regionflag == 0) {
        for (i = 0; i < nlocal; i++)
          if (mask[i] & groupbit) {
            if (rmass) values_one[bin[i]][m] += rmass[i];
            else values_one[bin[i]][m] += mass[type[i]];
          }
      } else {
        for (i = 0; i < nlocal; i++)
          if (mask[i] & groupbit && region->match(x[i][0],x[i][1],x[i][2])) {
            if (rmass) values_one[bin[i]][m] += rmass[i];
            else values_one[bin[i]][m] += mass[type[i]];
          }
      }

    // COMPUTE adds its scalar or vector component to values
    // invoke compute if not previously invoked

    } else if (which[m] == COMPUTE) {
      Compute *compute = modify->compute[n];
      if (!(compute->invoked_flag & INVOKED_PERATOM)) {
        compute->compute_peratom();
        compute->invoked_flag |= INVOKED_PERATOM;
      }
      double *vector = compute->vector_atom;
      double **array = compute->array_atom;
      int jm1 = j - 1;

      if (regionflag == 0) {
        for (i = 0; i < nlocal; i++)
          if (mask[i] & groupbit) {
            if (j == 0) values_one[bin[i]][m] += vector[i];
            else values_one[bin[i]][m] += array[i][jm1];
          }
      } else {
        for (i = 0; i < nlocal; i++)
          if (mask[i] & groupbit && region->match(x[i][0],x[i][1],x[i][2])) {
            if (j == 0) values_one[bin[i]][m] += vector[i];
            else values_one[bin[i]][m] += array[i][jm1];
          }
      }

    // FIX adds its scalar or vector component to values
    // access fix fields, guaranteed to be ready

    } else if (which[m] == FIX) {
      double *vector = modify->fix[n]->vector_atom;
      double **array = modify->fix[n]->array_atom;
      int jm1 = j - 1;

      if (regionflag == 0) {
        for (i = 0; i < nlocal; i++)
          if (mask[i] & groupbit) {
            if (j == 0) values_one[bin[i]][m] += vector[i];
            else values_one[bin[i]][m] += array[i][jm1];
          }
      } else {
        for (i = 0; i < nlocal; i++)
          if (mask[i] & groupbit && region->match(x[i][0],x[i][1],x[i][2])) {
            if (j == 0) values_one[bin[i]][m] += vector[i];
            else values_one[bin[i]][m] += array[i][jm1];
          }
      }

    // VARIABLE adds its per-atom quantities to values
    // evaluate atom-style variable

    } else if (which[m] == VARIABLE) {
      if (nlocal > maxvar) {
        maxvar = atom->nmax;
        memory->destroy(varatom);
        memory->create(varatom,maxvar,"ave/spatial:varatom");
      }

      input->variable->compute_atom(n,igroup,varatom,1,0);

      if (regionflag == 0) {
        for (i = 0; i < nlocal; i++)
          if (mask[i] & groupbit)
            values_one[bin[i]][m] += varatom[i];
      } else {
        for (i = 0; i < nlocal; i++)
          if (mask[i] & groupbit && region->match(x[i][0],x[i][1],x[i][2]))
            values_one[bin[i]][m] += varatom[i];
      }
    }
  }

  // process a single sample
  // if normflag = ALL, accumulate values,count separately to many
  // if normflag = SAMPLE, one = value/count, accumulate one to many
  // exception is SAMPLE density: no normalization by atom count

  if (normflag == ALL) {
    for (m = 0; m < nbins; m++) {
      count_many[m] += count_one[m];
      for (j = 0; j < nvalues; j++)
        values_many[m][j] += values_one[m][j];
    }
  } else {
    MPI_Allreduce(count_one,count_many,nbins,MPI_DOUBLE,MPI_SUM,world);
    for (m = 0; m < nbins; m++) {
      if (count_many[m] > 0.0)
        for (j = 0; j < nvalues; j++) {
          if (which[j] == DENSITY_NUMBER || which[j] == DENSITY_MASS)
            values_many[m][j] += values_one[m][j];
          else values_many[m][j] += values_one[m][j]/count_many[m];
        }
      count_sum[m] += count_many[m];
    }
  }

  // done if irepeat < nrepeat
  // else reset irepeat and nvalid

  irepeat++;
  if (irepeat < nrepeat) {
    nvalid += nevery;
    modify->addstep_compute(nvalid);
    return;
  }

  irepeat = 0;
  nvalid = ntimestep+nfreq - (nrepeat-1)*nevery;
  modify->addstep_compute(nvalid);

  // time average across samples
  // if normflag = ALL, final is total value / total count
  // if normflag = SAMPLE, final is sum of ave / repeat
  // exception is densities: normalized by repeat, not total count

  double repeat = nrepeat;
  double mv2d = force->mv2d;

  if (normflag == ALL)
  {
    MPI_Allreduce(count_many,count_sum,nbins,MPI_DOUBLE,MPI_SUM,world);
    MPI_Allreduce(&values_many[0][0],&values_sum[0][0],nbins*nvalues,
                  MPI_DOUBLE,MPI_SUM,world);
    for (m = 0; m < nbins; m++)
    {
      if (count_sum[m] > 0.0)
      {
        for (j = 0; j < nvalues; j++)
        {
          if (which[j] == DENSITY_NUMBER)
            values_sum[m][j] /= repeat;
          else if (which[j] == DENSITY_MASS)
            values_sum[m][j] *= mv2d/repeat;
          else
            values_sum[m][j] /= count_sum[m];
        }
      }
      count_sum[m] /= repeat;
    }
  }
  else
  {
    MPI_Allreduce(&values_many[0][0],&values_sum[0][0],nbins*nvalues,
                  MPI_DOUBLE,MPI_SUM,world);
    for (m = 0; m < nbins; m++)
    {
      for (j = 0; j < nvalues; j++)
        values_sum[m][j] /= repeat;
      count_sum[m] /= repeat;
    }
  }

  // density is additionally normalized by bin volume

  for (j = 0; j < nvalues; j++)
    if (which[j] == DENSITY_NUMBER || which[j] == DENSITY_MASS)
      for (m = 0; m < nbins; m++)
        values_sum[m][j] /= bin_volume;

  // if ave = ONE, only single Nfreq timestep value is needed
  // if ave = RUNNING, combine with all previous Nfreq timestep values
  // if ave = WINDOW, comine with nwindow most recent Nfreq timestep values

  if (ave == ONE) {
    for (m = 0; m < nbins; m++) {
      for (i = 0; i < nvalues; i++)
        values_total[m][i] = values_sum[m][i];
      count_total[m] = count_sum[m];
    }
    norm = 1;

  } else if (ave == RUNNING) {
    for (m = 0; m < nbins; m++) {
      for (i = 0; i < nvalues; i++)
        values_total[m][i] += values_sum[m][i];
      count_total[m] += count_sum[m];
    }
    norm++;

  } else if (ave == WINDOW) {
    for (m = 0; m < nbins; m++) {
      for (i = 0; i < nvalues; i++) {
        values_total[m][i] += values_sum[m][i];
        if (window_limit) values_total[m][i] -= values_list[iwindow][m][i];
        values_list[iwindow][m][i] = values_sum[m][i];
      }
      count_total[m] += count_sum[m];
      if (window_limit) count_total[m] -= count_list[iwindow][m];
      count_list[iwindow][m] = count_sum[m];
    }

    iwindow++;
    if (iwindow == nwindow) {
      iwindow = 0;
      window_limit = 1;
    }
    if (window_limit) norm = nwindow;
    else norm = iwindow;
  }

  // output result to file

  if (fp && me == 0) {
    if (overwrite) fseek(fp,filepos,SEEK_SET);
    if(write_ts_) fprintf(fp,BIGINT_FORMAT " %d\n",ntimestep,nbins);
    if (ndim == 1)
      for (m = 0; m < nbins; m++) {
        fprintf(fp,"  %d %g %g",m+1,coord[m][0],
                count_total[m]/norm);
        for (i = 0; i < nvalues; i++)
          fprintf(fp," %g",values_total[m][i]/norm);
        fprintf(fp,"\n");
      }
    else if (ndim == 2)
      for (m = 0; m < nbins; m++) {
        fprintf(fp,"  %d %g %g %g",m+1,coord[m][0],coord[m][1],
                count_total[m]/norm);
        for (i = 0; i < nvalues; i++)
          fprintf(fp," %g",values_total[m][i]/norm);
        fprintf(fp,"\n");
      }
    else
      for (m = 0; m < nbins; m++) {
        fprintf(fp,"  %d %g %g %g %g",m+1,coord[m][0],coord[m][1],coord[m][2],
                count_total[m]/norm);
        for (i = 0; i < nvalues; i++)
          fprintf(fp," %g",values_total[m][i]/norm);
        fprintf(fp,"\n");
      }

    fflush(fp);
    if (overwrite) {
      long fileend = ftell(fp);
      ftruncate(fileno(fp),fileend);
    }
  }

  // calculate mean and standard deviation
  
  if (calcStd != 0) {
    int n_samples = 0;
    int n_empty = 0;
    int n_lower = 0;
    int n_higher = 0;
    double *true_mean, *samples_mean, *std;
    double diff, count_sum, count_samples, samples_average, count_max;

    memory->create(true_mean,nvalues,"ave/spatial:true_mean");
    memory->create(samples_mean,nvalues,"ave/spatial:samples_mean");
    memory->create(std,nvalues,"ave/spatial:std");

    count_sum = 0;
    count_samples = 0;
    samples_average = 0;
    count_max = 0;

    for (i = 0; i < nvalues; i++) {
      true_mean[i] = 0;

      for (m = 0; m < nbins; m++) {
        true_mean[i] += values_total[m][i] * count_total[m];
        if (i == 0) count_sum += count_total[m];
      }

      if (count_sum != 0) true_mean[i] /= count_sum;

      samples_mean[i] = 0;
      std[i] = 0;
    }

    for (m = 0; m < nbins; m++) {
      if (count_total[m] > count_max) count_max = count_total[m];
      if (count_total[m] == 0) n_empty += 1;
      else if (count_total[m] < lowerLimit) n_lower += 1;
      else if (count_total[m] > upperLimit) n_higher += 1;
      else {
        n_samples += 1;
        count_samples += count_total[m];
        for (i = 0; i < nvalues; i++) {
          samples_mean[i] += values_total[m][i];
          diff = values_total[m][i] - true_mean[i];
          std[i] += diff * diff;
        }
      }
    }

    if (n_samples != 0) {
      for (i = 0; i < nvalues; i++) {
        samples_average = count_samples / n_samples; // average # particles per sample
        samples_mean[i] /= n_samples;
        std[i] /= n_samples; // not(n_sample-1) since true_mean is the true Average
        std[i] = sqrt(std[i]);
        samples_mean[i] /= norm;
        std[i] /= norm;
      }
    }

    // output to file
    if (fp2 && me == 0) {
      fprintf(fp2,BIGINT_FORMAT " %g %d %g %d %d %d %d %g",
        ntimestep,count_sum,nbins,count_max,n_empty,n_lower,n_higher,n_samples,samples_average);
      for (i = 0; i < nvalues; i++) {
        fprintf(fp2," %g %g %g",true_mean[i],samples_mean[i],std[i]); //;
      }
      fprintf(fp2,"\n");
      fflush(fp2);
    }

    memory->destroy(true_mean);
    memory->destroy(samples_mean);
    memory->destroy(std);
  }
}

/* ----------------------------------------------------------------------
   setup 1d, 2d, or 3d bins and their extent and coordinates
   called at setup() and when averaging occurs if box size changes
------------------------------------------------------------------------- */

void FixAveSpatial::setup_bins()
{
  int i,j,k,m,n;
  double lo,hi,coord1,coord2;

  // lo = bin boundary immediately below boxlo
  // hi = bin boundary immediately above boxhi
  // allocate and initialize arrays based on new bin count

  double *boxlo,*boxhi,*prd;
  if (scaleflag == REDUCED) {
    boxlo = domain->boxlo_lamda;
    boxhi = domain->boxhi_lamda;
    prd = domain->prd_lamda;
  } else {
    boxlo = domain->boxlo;
    boxhi = domain->boxhi;
    prd = domain->prd;
  }

  if (domain->dimension == 3)
    bin_volume = domain->xprd * domain->yprd * domain->zprd;
  else bin_volume = domain->xprd * domain->yprd;
  nbins = 1;

  for (m = 0; m < ndim; m++) {
    if (originflag[m] == LOWER) origin[m] = boxlo[dim[m]];
    else if (originflag[m] == UPPER) origin[m] = boxhi[dim[m]];
    else if (originflag[m] == CENTER)
      origin[m] = 0.5 * (boxlo[dim[m]] + boxhi[dim[m]]);

    if (origin[m] < boxlo[dim[m]]) {
      n = static_cast<int> ((boxlo[dim[m]] - origin[m]) * invdelta[m]);
      lo = origin[m] + n*delta[m];
    } else {
      n = static_cast<int> ((origin[m] - boxlo[dim[m]]) * invdelta[m]);
      lo = origin[m] - n*delta[m];
      if (lo > boxlo[dim[m]]) lo -= delta[m];
    }
    if (origin[m] < boxhi[dim[m]]) {
      n = static_cast<int> ((boxhi[dim[m]] - origin[m]) * invdelta[m]);
      hi = origin[m] + n*delta[m];
      if (hi < boxhi[dim[m]]) hi += delta[m];
    } else {
      n = static_cast<int> ((origin[m] - boxhi[dim[m]]) * invdelta[m]);
      hi = origin[m] - n*delta[m];
    }

    offset[m] = lo;
    nlayers[m] = static_cast<int> ((hi-lo) * invdelta[m] + 0.5);
    nbins *= nlayers[m];
    bin_volume *= delta[m]/prd[dim[m]];
  }

  // reallocate bin arrays if needed

  if (nbins > maxbin) {
    maxbin = nbins;
    memory->grow(count_one,nbins,"ave/spatial:count_one");
    memory->grow(count_many,nbins,"ave/spatial:count_many");
    memory->grow(count_sum,nbins,"ave/spatial:count_sum");
    memory->grow(count_total,nbins,"ave/spatial:count_total");

    memory->grow(coord,nbins,ndim,"ave/spatial:coord");
    memory->grow(values_one,nbins,nvalues,"ave/spatial:values_one");
    memory->grow(values_many,nbins,nvalues,"ave/spatial:values_many");
    memory->grow(values_sum,nbins,nvalues,"ave/spatial:values_sum");
    memory->grow(values_total,nbins,nvalues,"ave/spatial:values_total");

    // only allocate count and values list for ave = WINDOW

    if (ave == WINDOW) {
      memory->create(count_list,nwindow,nbins,"ave/spatial:count_list");
      memory->create(values_list,nwindow,nbins,nvalues,
                     "ave/spatial:values_list");
    }

    // reinitialize regrown count/values total since they accumulate

    for (m = 0; m < nbins; m++) {
      for (i = 0; i < nvalues; i++) values_total[m][i] = 0.0;
      count_total[m] = 0.0;
    }
  }

  // set bin coordinates

  if (ndim == 1) {
    for (i = 0; i < nlayers[0]; i++)
      coord[i][0] = offset[0] + (i+0.5)*delta[0];
  } else if (ndim == 2) {
    m = 0;
    for (i = 0; i < nlayers[0]; i++) {
      coord1 = offset[0] + (i+0.5)*delta[0];
      for (j = 0; j < nlayers[1]; j++) {
        coord[m][0] = coord1;
        coord[m][1] = offset[1] + (j+0.5)*delta[1];
        m++;
      }
    }
  } else if (ndim == 3) {
    m = 0;
    for (i = 0; i < nlayers[0]; i++) {
      coord1 = offset[0] + (i+0.5)*delta[0];
      for (j = 0; j < nlayers[1]; j++) {
        coord2 = offset[1] + (j+0.5)*delta[1];
        for (k = 0; k < nlayers[2]; k++) {
          coord[m][0] = coord1;
          coord[m][1] = coord2;
          coord[m][2] = offset[2] + (k+0.5)*delta[2];
          m++;
        }
      }
    }
  }
}

/* ----------------------------------------------------------------------
   assign each atom to a 1d bin
------------------------------------------------------------------------- */

void FixAveSpatial::atom2bin1d()
{
  int i,ibin;
  double *boxlo=NULL,*boxhi=NULL,*prd=NULL;
  double xremap;
  double lamda[3];

  double **x = atom->x;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  int idim = dim[0];
  int nlayerm1 = nlayers[0] - 1;
  int periodicity = domain->periodicity[idim];

  if (periodicity) {
    if (scaleflag == REDUCED) {
      boxlo = domain->boxlo_lamda;
      boxhi = domain->boxhi_lamda;
      prd = domain->prd_lamda;
    } else {
      boxlo = domain->boxlo;
      boxhi = domain->boxhi;
      prd = domain->prd;
    }
  }

  // remap each atom's relevant coord back into box via PBC if necessary
  // if scaleflag = REDUCED, box coords -> lamda coords

  if (regionflag == 0) {
    if (scaleflag == REDUCED) domain->x2lamda(nlocal);
    for (i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        xremap = x[i][idim];
        if (periodicity) {
          if (xremap < boxlo[idim]) xremap += prd[idim];
          if (xremap >= boxhi[idim]) xremap -= prd[idim];
        }
        ibin = static_cast<int> ((xremap - offset[0]) * invdelta[0]);
        ibin = MAX(ibin,0);
        ibin = MIN(ibin,nlayerm1);
        bin[i] = ibin;
        count_one[ibin] += 1.0;
      }
    if (scaleflag == REDUCED) domain->lamda2x(nlocal);

  } else {
    for (i = 0; i < nlocal; i++)
      if (mask[i] & groupbit && region->match(x[i][0],x[i][1],x[i][2])) {
        if (scaleflag == REDUCED) {
          domain->x2lamda(x[i],lamda);
          xremap = lamda[idim];
        } else xremap = x[i][idim];
        if (periodicity) {
          if (xremap < boxlo[idim]) xremap += prd[idim];
          if (xremap >= boxhi[idim]) xremap -= prd[idim];
        }
        ibin = static_cast<int> ((xremap - offset[0]) * invdelta[0]);
        ibin = MAX(ibin,0);
        ibin = MIN(ibin,nlayerm1);
        bin[i] = ibin;
        count_one[ibin] += 1.0;
      }
  }
}

/* ----------------------------------------------------------------------
   assign each atom to a 2d bin
------------------------------------------------------------------------- */

void FixAveSpatial::atom2bin2d()
{
  int i,ibin,i1bin,i2bin;
  double *boxlo=NULL,*boxhi=NULL,*prd=NULL;
  double xremap,yremap;
  double lamda[3];

  double **x = atom->x;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  int idim = dim[0];
  int jdim = dim[1];
  int nlayer1m1 = nlayers[0] - 1;
  int nlayer2m1 = nlayers[1] - 1;
  int* periodicity = domain->periodicity;

  if (periodicity[idim] || periodicity[jdim]) {
    if (scaleflag == REDUCED) {
      boxlo = domain->boxlo_lamda;
      boxhi = domain->boxhi_lamda;
      prd = domain->prd_lamda;
    } else {
      boxlo = domain->boxlo;
      boxhi = domain->boxhi;
      prd = domain->prd;
    }
  }

  // remap each atom's relevant coord back into box via PBC if necessary
  // if scaleflag = REDUCED, box coords -> lamda coords

  if (regionflag == 0) {
    if (scaleflag == REDUCED) domain->x2lamda(nlocal);
    for (i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        xremap = x[i][idim];
        if (periodicity[idim]) {
          if (xremap < boxlo[idim]) xremap += prd[idim];
          if (xremap >= boxhi[idim]) xremap -= prd[idim];
        }
        i1bin = static_cast<int> ((xremap - offset[0]) * invdelta[0]);
        i1bin = MAX(i1bin,0);
        i1bin = MIN(i1bin,nlayer1m1);

        yremap = x[i][jdim];
        if (periodicity[jdim]) {
          if (yremap < boxlo[jdim]) yremap += prd[jdim];
          if (yremap >= boxhi[jdim]) yremap -= prd[jdim];
        }
        i2bin = static_cast<int> ((yremap - offset[1]) * invdelta[1]);
        i2bin = MAX(i2bin,0);
        i2bin = MIN(i2bin,nlayer2m1);

        ibin = i1bin*nlayers[1] + i2bin;
        bin[i] = ibin;
        count_one[ibin] += 1.0;
      }
    if (scaleflag == REDUCED) domain->lamda2x(nlocal);

  } else {
    for (i = 0; i < nlocal; i++)
      if (mask[i] & groupbit && region->match(x[i][0],x[i][1],x[i][2])) {
        if (scaleflag == REDUCED) {
          domain->x2lamda(x[i],lamda);
          xremap = lamda[idim];
          yremap = lamda[jdim];
        } else {
          xremap = x[i][idim];
          yremap = x[i][jdim];
        }

        if (periodicity[idim]) {
          if (xremap < boxlo[idim]) xremap += prd[idim];
          if (xremap >= boxhi[idim]) xremap -= prd[idim];
        }
        i1bin = static_cast<int> ((xremap - offset[0]) * invdelta[0]);
        i1bin = MAX(i1bin,0);
        i1bin = MIN(i1bin,nlayer1m1-1);

        if (periodicity[jdim]) {
          if (yremap < boxlo[jdim]) yremap += prd[jdim];
          if (yremap >= boxhi[jdim]) yremap -= prd[jdim];
        }
        i2bin = static_cast<int> ((yremap - offset[1]) * invdelta[1]);
        i2bin = MAX(i2bin,0);
        i2bin = MIN(i2bin,nlayer2m1-1);

        ibin = i1bin*nlayers[1] + i2bin;
        bin[i] = ibin;
        count_one[ibin] += 1.0;
      }
  }
}

/* ----------------------------------------------------------------------
   assign each atom to a 3d bin
------------------------------------------------------------------------- */

void FixAveSpatial::atom2bin3d()
{
  int i,ibin,i1bin,i2bin,i3bin;
  double *boxlo=NULL,*boxhi=NULL,*prd=NULL;
  double xremap,yremap,zremap;
  double lamda[3];

  double **x = atom->x;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  int idim = dim[0];
  int jdim = dim[1];
  int kdim = dim[2];
  int nlayer1m1 = nlayers[0] - 1;
  int nlayer2m1 = nlayers[1] - 1;
  int nlayer3m1 = nlayers[2] - 1;
  int* periodicity = domain->periodicity;

  if (periodicity[idim] || periodicity[jdim] || periodicity[kdim]) {
    if (scaleflag == REDUCED) {
      boxlo = domain->boxlo_lamda;
      boxhi = domain->boxhi_lamda;
      prd = domain->prd_lamda;
    } else {
      boxlo = domain->boxlo;
      boxhi = domain->boxhi;
      prd = domain->prd;
    }
  }

  // remap each atom's relevant coord back into box via PBC if necessary
  // if scaleflag = REDUCED, box coords -> lamda coords

  if (regionflag == 0) {
    if (scaleflag == REDUCED) domain->x2lamda(nlocal);
    for (i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        xremap = x[i][idim];
        if (periodicity[idim]) {
          if (xremap < boxlo[idim]) xremap += prd[idim];
          if (xremap >= boxhi[idim]) xremap -= prd[idim];
        }
        i1bin = static_cast<int> ((xremap - offset[0]) * invdelta[0]);
        i1bin = MAX(i1bin,0);
        i1bin = MIN(i1bin,nlayer1m1);

        yremap = x[i][jdim];
        if (periodicity[jdim]) {
          if (yremap < boxlo[jdim]) yremap += prd[jdim];
          if (yremap >= boxhi[jdim]) yremap -= prd[jdim];
        }
        i2bin = static_cast<int> ((yremap - offset[1]) * invdelta[1]);
        i2bin = MAX(i2bin,0);
        i2bin = MIN(i2bin,nlayer2m1);

        zremap = x[i][kdim];
        if (periodicity[kdim]) {
          if (zremap < boxlo[kdim]) zremap += prd[kdim];
          if (zremap >= boxhi[kdim]) zremap -= prd[kdim];
        }
        i3bin = static_cast<int> ((zremap - offset[2]) * invdelta[2]);
        i3bin = MAX(i3bin,0);
        i3bin = MIN(i3bin,nlayer3m1);

        ibin = i1bin*nlayers[1]*nlayers[2] + i2bin*nlayers[2] + i3bin;
        bin[i] = ibin;
        count_one[ibin] += 1.0;
      }
    if (scaleflag == REDUCED) domain->lamda2x(nlocal);

  } else {
    for (i = 0; i < nlocal; i++)
      if (mask[i] & groupbit && region->match(x[i][0],x[i][1],x[i][2])) {
        if (scaleflag == REDUCED) {
          domain->x2lamda(x[i],lamda);
          xremap = lamda[idim];
          yremap = lamda[jdim];
          zremap = lamda[kdim];
        } else {
          xremap = x[i][idim];
          yremap = x[i][jdim];
          zremap = x[i][kdim];
        }

        if (periodicity[idim]) {
          if (xremap < boxlo[idim]) xremap += prd[idim];
          if (xremap >= boxhi[idim]) xremap -= prd[idim];
        }
        i1bin = static_cast<int> ((xremap - offset[0]) * invdelta[0]);
        i1bin = MAX(i1bin,0);
        i1bin = MIN(i1bin,nlayer1m1);

        if (periodicity[jdim]) {
          if (yremap < boxlo[jdim]) yremap += prd[jdim];
          if (yremap >= boxhi[jdim]) yremap -= prd[jdim];
        }
        i2bin = static_cast<int> ((yremap - offset[1]) * invdelta[1]);
        i2bin = MAX(i2bin,0);
        i2bin = MIN(i2bin,nlayer2m1);

        if (periodicity[kdim]) {
          if (zremap < boxlo[kdim]) zremap += prd[kdim];
          if (zremap >= boxhi[kdim]) zremap -= prd[kdim];
        }
        i3bin = static_cast<int> ((zremap - offset[2]) * invdelta[2]);
        i3bin = MAX(i3bin,0);
        i3bin = MIN(i3bin,nlayer3m1);

        ibin = i1bin*nlayers[1]*nlayers[2] + i2bin*nlayers[2] + i3bin;
        bin[i] = ibin;
        count_one[ibin] += 1.0;
      }
  }
}

/* ----------------------------------------------------------------------
   return I,J array value
   if I exceeds current bins, return 0.0 instead of generating an error
   column 1,2,3 = bin coords, next column = count, remaining columns = Nvalues
------------------------------------------------------------------------- */

double FixAveSpatial::compute_array(int i, int j)
{
  if (values_total == NULL) return 0.0;
  if (i >= nbins) return 0.0;
  if (j < ndim) return coord[i][j];
  j -= ndim+1;
  if (!norm) return 0.0;
  if (j < 0) return count_total[i]/norm;
  return values_total[i][j]/norm;
}

/* ----------------------------------------------------------------------
   calculate nvalid = next step on which end_of_step does something
   can be this timestep if multiple of nfreq and nrepeat = 1
   else backup from next multiple of nfreq
------------------------------------------------------------------------- */

bigint FixAveSpatial::nextvalid()
{
  bigint nvalid = (update->ntimestep/nfreq)*nfreq + nfreq;
  if (nvalid-nfreq == update->ntimestep && nrepeat == 1)
    nvalid = update->ntimestep;
  else
    nvalid -= (nrepeat-1)*nevery;
  if (nvalid < update->ntimestep) nvalid += nfreq;
  return nvalid;
}

/* ----------------------------------------------------------------------
   memory usage of varatom and bins
------------------------------------------------------------------------- */

double FixAveSpatial::memory_usage()
{
  double bytes = maxvar * sizeof(double);         // varatom
  bytes += maxatom * sizeof(int);                 // bin
  bytes += 4*nbins * sizeof(double);              // count one,many,sum,total
  bytes += ndim*nbins * sizeof(double);           // coord
  bytes += nvalues*nbins * sizeof(double);        // values one,many,sum,total
  bytes += nwindow*nbins * sizeof(double);          // count_list
  bytes += nwindow*nbins*nvalues * sizeof(double);  // values_list
  return bytes;
}

/* ---------------------------------------------------------------------- */

void FixAveSpatial::reset_timestep(bigint ntimestep)
{
  if (ntimestep > nvalid) error->all(FLERR,"Fix ave/spatial missed timestep");
}
