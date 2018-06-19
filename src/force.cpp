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

#include <stdlib.h>
#include <string.h>
#include "ctype.h"
#include "force.h"
#include "style_bond.h"
#include "style_angle.h"
#include "style_dihedral.h"
#include "style_improper.h"
#include "style_pair.h"
#include "style_kspace.h"
#include "atom.h"
#include "comm.h"
#include "pair.h"
#include "pair_gran.h"
#include "pair_hybrid.h"
#include "pair_hybrid_overlay.h"
#include "bond.h"
#include "bond_hybrid.h"
#include "angle.h"
#include "dihedral.h"
#include "improper.h"
#include "kspace.h"
#include "group.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

Force::Force(LAMMPS *lmp) : Pointers(lmp), registry(lmp)
{
  newton = newton_pair = newton_bond = 1;

  special_lj[0] = special_coul[0] = 1.0;
  special_lj[1] = special_lj[2] = special_lj[3] = 0.0;
  special_coul[1] = special_coul[2] = special_coul[3] = 0.0;
  special_angle = special_dihedral = 0;
  special_extra = 0;

  dielectric = 1.0;

  coarsegraining_ = 1.0;
  error_coarsegraining_ = true;
  warn_coarsegraining_ = false;

  pair = NULL;
  bond = NULL;
  angle = NULL;
  dihedral = NULL;
  improper = NULL;
  kspace = NULL;

  char *str = (char *) "none";
  int n = strlen(str) + 1;
  pair_style = new char[n];
  strcpy(pair_style,str);
  bond_style = new char[n];
  strcpy(bond_style,str);
  angle_style = new char[n];
  strcpy(angle_style,str);
  dihedral_style = new char[n];
  strcpy(dihedral_style,str);
  improper_style = new char[n];
  strcpy(improper_style,str);
  kspace_style = new char[n];
  strcpy(kspace_style,str);

  // fill pair map with pair styles listed in style_pair.h

  pair_map = new std::map<std::string,PairCreator>();

#define PAIR_CLASS
#define PairStyle(key,Class) \
  (*pair_map)[#key] = &pair_creator<Class>;
#include "style_pair.h"
#undef PairStyle
#undef PAIR_CLASS
}

/* ---------------------------------------------------------------------- */

Force::~Force()
{
  delete [] pair_style;
  delete [] bond_style;
  delete [] angle_style;
  delete [] dihedral_style;
  delete [] improper_style;
  delete [] kspace_style;

  if (pair) delete pair;
  if (bond) delete bond;
  if (angle) delete angle;
  if (dihedral) delete dihedral;
  if (improper) delete improper;
  if (kspace) delete kspace;

  delete pair_map;
}

/* ---------------------------------------------------------------------- */

void Force::init()
{
  registry.init();
  qqrd2e = qqr2e/dielectric;

  if (kspace) kspace->init();         // kspace must come before pair
  if (pair) pair->init();             // so g_ewald is defined
  if (bond) bond->init();
  if (angle) angle->init();
  if (dihedral) dihedral->init();
  if (improper) improper->init();
  if(cg_active() && warn_cg() && atom->ntypes != int(coarsegrainingTypeBased_.size()))
    error->warningAll(FLERR,"Coarse graining factor not specified for all atom types. will use maximum CG for unspecified atom types.\n\n");
}

/* ----------------------------------------------------------------------
   create a pair style, called from input script
------------------------------------------------------------------------- */

void Force::create_pair(const char *style, const char *suffix)
{
  delete [] pair_style;
  if (pair) delete pair;

  int sflag;
  pair = new_pair(style,suffix,sflag);

  if (sflag) {
    char estyle[256];
    sprintf(estyle,"%s/%s",style,suffix);
    int n = strlen(estyle) + 1;
    pair_style = new char[n];
    strcpy(pair_style,estyle);
  } else {
    int n = strlen(style) + 1;
    pair_style = new char[n];
    strcpy(pair_style,style);
  }
}

/* ----------------------------------------------------------------------
   create a pair style, called from restart file
------------------------------------------------------------------------- */

void Force::create_pair_from_restart(FILE * fp, const char *style, const char *suffix)
{
  delete [] pair_style;
  if (pair) delete pair;

  int sflag;
  pair = new_pair_from_restart(fp, style,suffix,sflag);

  if (sflag) {
    char estyle[256];
    sprintf(estyle,"%s/%s",style,suffix);
    int n = strlen(estyle) + 1;
    pair_style = new char[n];
    strcpy(pair_style,estyle);
  } else {
    int n = strlen(style) + 1;
    pair_style = new char[n];
    strcpy(pair_style,style);
  }
}

/* ----------------------------------------------------------------------
   generate a pair class
   try first with suffix appended
------------------------------------------------------------------------- */

Pair *Force::new_pair(const char *style, const char *suffix, int &sflag)
{
  if (suffix && lmp->suffix_enable) {
    sflag = 1;
    char estyle[256];
    sprintf(estyle,"%s/%s",style,suffix);

    if (pair_map->find(estyle) != pair_map->end()) {
      PairCreator pair_creator = (*pair_map)[estyle];
      return pair_creator(lmp);
    }

  }

  sflag = 0;

  if (strcmp(style,"none") == 0) return NULL;
  if (pair_map->find(style) != pair_map->end()) {
    PairCreator pair_creator = (*pair_map)[style];
    return pair_creator(lmp);
  }

  error->all(FLERR,"Invalid pair style");

  return NULL;
}

/* ----------------------------------------------------------------------
   generate a pair class from restart file, first with suffix appended
------------------------------------------------------------------------- */

Pair *Force::new_pair_from_restart(FILE * fp, const char *style, const char *suffix, int &sflag)
{
  if (suffix && lmp->suffix_enable) {
    sflag = 1;
    char estyle[256];
    sprintf(estyle,"%s/%s",style,suffix);

    if (pair_map->find(estyle) != pair_map->end()) {
      PairCreator pair_creator = (*pair_map)[estyle];
      return pair_creator(lmp);
    }
  }

  sflag = 0;

  if (strcmp(style,"none") == 0) return NULL;
  if (pair_map->find(style) != pair_map->end()) {
    PairCreator pair_creator = (*pair_map)[style];
    return pair_creator(lmp);
  }

  error->all(FLERR,"Invalid pair style");
  return NULL;
}

/* ----------------------------------------------------------------------
   one instance per pair style in style_pair.h
------------------------------------------------------------------------- */

template <typename T>
Pair *Force::pair_creator(LAMMPS *lmp)
{
  return new T(lmp);
}

/* ----------------------------------------------------------------------
   return ptr to Pair class if matches word or matches hybrid sub-style
   if exact, then style name must be exact match to word
   if not exact, style name must contain word
   return NULL if no match or multiple sub-styles match
------------------------------------------------------------------------- */

Pair *Force::pair_match(const char *word, int exact)
{
  int iwhich,count;

  if (exact && strcmp(pair_style,word) == 0) return pair;
  else if (!exact && strstr(pair_style,word)) return pair;

  else if (strstr(pair_style,"hybrid/overlay")) {
    PairHybridOverlay *hybrid = (PairHybridOverlay *) pair;
    count = 0;
    for (int i = 0; i < hybrid->nstyles; i++)
      if ((exact && strcmp(hybrid->keywords[i],word) == 0) ||
          (!exact && strstr(hybrid->keywords[i],word))) {
        iwhich = i;
        count++;
      }
    if (count == 1) return hybrid->styles[iwhich];

  } else if (strstr(pair_style,"hybrid")) {
    PairHybrid *hybrid = (PairHybrid *) pair;
    count = 0;
    for (int i = 0; i < hybrid->nstyles; i++)
      if ((exact && strcmp(hybrid->keywords[i],word) == 0) ||
          (!exact && strstr(hybrid->keywords[i],word))) {
        iwhich = i;
        count++;
      }
    if (count == 1) return hybrid->styles[iwhich];
  }

  return NULL;
}

/* ----------------------------------------------------------------------
   create a bond style, called from input script or restart file
------------------------------------------------------------------------- */

void Force::create_bond(const char *style, const char *suffix)
{
  delete [] bond_style;
  if (bond) delete bond;

  int sflag;
  bond = new_bond(style,suffix,sflag);

  if (sflag) {
    char estyle[256];
    sprintf(estyle,"%s/%s",style,suffix);
    int n = strlen(estyle) + 1;
    bond_style = new char[n];
    strcpy(bond_style,estyle);
  } else {
    int n = strlen(style) + 1;
    bond_style = new char[n];
    strcpy(bond_style,style);
  }
}

/* ----------------------------------------------------------------------
   generate a bond class, fist with suffix appended
------------------------------------------------------------------------- */

Bond *Force::new_bond(const char *style, const char *suffix, int &sflag)
{
  if (suffix && lmp->suffix_enable) {
    sflag = 1;
    char estyle[256];
    sprintf(estyle,"%s/%s",style,suffix);

    if (0) return NULL;

#define BOND_CLASS
#define BondStyle(key,Class) \
    else if (strcmp(estyle,#key) == 0) return new Class(lmp);
#include "style_bond.h"
#undef BondStyle
#undef BOND_CLASS
  }

  sflag = 0;

  if (strcmp(style,"none") == 0) return NULL;

#define BOND_CLASS
#define BondStyle(key,Class) \
  else if (strcmp(style,#key) == 0) return new Class(lmp);
#include "style_bond.h"
#undef BOND_CLASS

  else error->all(FLERR,"Invalid bond style");
  return NULL;
}

/* ----------------------------------------------------------------------
   return ptr to current bond class or hybrid sub-class if matches style
------------------------------------------------------------------------- */

Bond *Force::bond_match(const char *style)
{
  if (strcmp(bond_style,style) == 0) return bond;
  else if (strcmp(bond_style,"hybrid") == 0) {
    BondHybrid *hybrid = (BondHybrid *) bond;
    for (int i = 0; i < hybrid->nstyles; i++)
      if (strcmp(hybrid->keywords[i],style) == 0) return hybrid->styles[i];
  }
  return NULL;
}

/* ----------------------------------------------------------------------
   create an angle style, called from input script or restart file
------------------------------------------------------------------------- */

void Force::create_angle(const char *style, const char *suffix)
{
  delete [] angle_style;
  if (angle) delete angle;

  int sflag;
  angle = new_angle(style,suffix,sflag);

  if (sflag) {
    char estyle[256];
    sprintf(estyle,"%s/%s",style,suffix);
    int n = strlen(estyle) + 1;
    angle_style = new char[n];
    strcpy(angle_style,estyle);
  } else {
    int n = strlen(style) + 1;
    angle_style = new char[n];
    strcpy(angle_style,style);
  }
}

/* ----------------------------------------------------------------------
   generate an angle class
------------------------------------------------------------------------- */

Angle *Force::new_angle(const char *style, const char *suffix, int &sflag)
{
  if (suffix && lmp->suffix_enable) {
    sflag = 1;
    char estyle[256];
    sprintf(estyle,"%s/%s",style,suffix);

    if (0) return NULL;

#define ANGLE_CLASS
#define AngleStyle(key,Class) \
    else if (strcmp(estyle,#key) == 0) return new Class(lmp);
#include "style_angle.h"
#undef AngleStyle
#undef ANGLE_CLASS

  }

  sflag = 0;

  if (strcmp(style,"none") == 0) return NULL;

#define ANGLE_CLASS
#define AngleStyle(key,Class) \
  else if (strcmp(style,#key) == 0) return new Class(lmp);
#include "style_angle.h"
#undef ANGLE_CLASS

  else error->all(FLERR,"Invalid angle style");
  return NULL;
}

/* ----------------------------------------------------------------------
   create a dihedral style, called from input script or restart file
------------------------------------------------------------------------- */

void Force::create_dihedral(const char *style, const char *suffix)
{
  delete [] dihedral_style;
  if (dihedral) delete dihedral;

  int sflag;
  dihedral = new_dihedral(style,suffix,sflag);

  if (sflag) {
    char estyle[256];
    sprintf(estyle,"%s/%s",style,suffix);
    int n = strlen(estyle) + 1;
    dihedral_style = new char[n];
    strcpy(dihedral_style,estyle);
  } else {
    int n = strlen(style) + 1;
    dihedral_style = new char[n];
    strcpy(dihedral_style,style);
  }
}

/* ----------------------------------------------------------------------
   generate a dihedral class
------------------------------------------------------------------------- */

Dihedral *Force::new_dihedral(const char *style, const char *suffix, int &sflag)
{
  if (suffix && lmp->suffix_enable) {
    sflag = 1;
    char estyle[256];
    sprintf(estyle,"%s/%s",style,suffix);

    if (0) return NULL;

#define DIHEDRAL_CLASS
#define DihedralStyle(key,Class) \
    else if (strcmp(estyle,#key) == 0) return new Class(lmp);
#include "style_dihedral.h"
#undef DihedralStyle
#undef DIHEDRAL_CLASS

  }

  sflag = 0;

  if (strcmp(style,"none") == 0) return NULL;

#define DIHEDRAL_CLASS
#define DihedralStyle(key,Class) \
  else if (strcmp(style,#key) == 0) return new Class(lmp);
#include "style_dihedral.h"
#undef DihedralStyle
#undef DIHEDRAL_CLASS

  else error->all(FLERR,"Invalid dihedral style");
  return NULL;
}

/* ----------------------------------------------------------------------
   create an improper style, called from input script or restart file
------------------------------------------------------------------------- */

void Force::create_improper(const char *style, const char *suffix)
{
  delete [] improper_style;
  if (improper) delete improper;

  int sflag;
  improper = new_improper(style,suffix,sflag);

  if (sflag) {
    char estyle[256];
    sprintf(estyle,"%s/%s",style,suffix);
    int n = strlen(estyle) + 1;
    improper_style = new char[n];
    strcpy(improper_style,estyle);
  } else {
    int n = strlen(style) + 1;
    improper_style = new char[n];
    strcpy(improper_style,style);
  }
}

/* ----------------------------------------------------------------------
   generate a improper class
------------------------------------------------------------------------- */

Improper *Force::new_improper(const char *style, const char *suffix, int &sflag)
{
  if (suffix && lmp->suffix_enable) {
    sflag = 1;
    char estyle[256];
    sprintf(estyle,"%s/%s",style,suffix);

    if (0) return NULL;

#define IMPROPER_CLASS
#define ImproperStyle(key,Class) \
    else if (strcmp(estyle,#key) == 0) return new Class(lmp);
#include "style_improper.h"
#undef ImproperStyle
#undef IMPROPER_CLASS

  }

  sflag = 0;

  if (strcmp(style,"none") == 0) return NULL;

#define IMPROPER_CLASS
#define ImproperStyle(key,Class) \
  else if (strcmp(style,#key) == 0) return new Class(lmp);
#include "style_improper.h"
#undef IMPROPER_CLASS

  else error->all(FLERR,"Invalid improper style");
  return NULL;
}

/* ----------------------------------------------------------------------
   new kspace style
------------------------------------------------------------------------- */

void Force::create_kspace(int narg, char **arg, const char *suffix)
{
  delete [] kspace_style;
  if (kspace) delete kspace;

  int sflag;
  kspace = new_kspace(narg,arg,suffix,sflag);

  if (sflag) {
    char estyle[256];
    sprintf(estyle,"%s/%s",arg[0],suffix);
    int n = strlen(estyle) + 1;
    kspace_style = new char[n];
    strcpy(kspace_style,estyle);
  } else {
    int n = strlen(arg[0]) + 1;
    kspace_style = new char[n];
    strcpy(kspace_style,arg[0]);
  }
}

/* ----------------------------------------------------------------------
   generate a kspace class
------------------------------------------------------------------------- */

KSpace *Force::new_kspace(int narg, char **arg, const char *suffix, int &sflag)
{
  if (suffix && lmp->suffix_enable) {
    sflag = 1;
    char estyle[256];
    sprintf(estyle,"%s/%s",arg[0],suffix);

    if (0) return NULL;

#define KSPACE_CLASS
#define KSpaceStyle(key,Class) \
  else if (strcmp(estyle,#key) == 0) return new Class(lmp,narg-1,&arg[1]);
#include "style_kspace.h"
#undef KSpaceStyle
#undef KSPACE_CLASS

  }

  sflag = 0;

  if (strcmp(arg[0],"none") == 0) return NULL;

#define KSPACE_CLASS
#define KSpaceStyle(key,Class) \
  else if (strcmp(arg[0],#key) == 0) return  new Class(lmp,narg-1,&arg[1]);
#include "style_kspace.h"
#undef KSPACE_CLASS

  else error->all(FLERR,"Invalid kspace style");
  return NULL;
}

/* ----------------------------------------------------------------------
   return ptr to Kspace class if matches word
   if exact, then style name must be exact match to word
   if not exact, style name must contain word
   return NULL if no match
------------------------------------------------------------------------- */

KSpace *Force::kspace_match(const char *word, int exact)
{
  if (exact && strcmp(kspace_style,word) == 0) return kspace;
  else if (!exact && strstr(kspace_style,word)) return kspace;
  return NULL;
}

/* ----------------------------------------------------------------------
   set special bond values
------------------------------------------------------------------------- */

void Force::set_special(int narg, char **arg)
{
  if (narg == 0) error->all(FLERR,"Illegal special_bonds command");

  int iarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"amber") == 0) {
      if (iarg+1 > narg) error->all(FLERR,"Illegal special_bonds command");
      special_lj[1] = 0.0;
      special_lj[2] = 0.0;
      special_lj[3] = 0.5;
      special_coul[1] = 0.0;
      special_coul[2] = 0.0;
      special_coul[3] = 5.0/6.0;
      iarg += 1;
    } else if (strcmp(arg[iarg],"charmm") == 0) {
      if (iarg+1 > narg) error->all(FLERR,"Illegal special_bonds command");
      special_lj[1] = 0.0;
      special_lj[2] = 0.0;
      special_lj[3] = 0.0;
      special_coul[1] = 0.0;
      special_coul[2] = 0.0;
      special_coul[3] = 0.0;
      iarg += 1;
    } else if (strcmp(arg[iarg],"dreiding") == 0) {
      if (iarg+1 > narg) error->all(FLERR,"Illegal special_bonds command");
      special_lj[1] = 0.0;
      special_lj[2] = 0.0;
      special_lj[3] = 1.0;
      special_coul[1] = 0.0;
      special_coul[2] = 0.0;
      special_coul[3] = 1.0;
      iarg += 1;
    } else if (strcmp(arg[iarg],"fene") == 0) {
      if (iarg+1 > narg) error->all(FLERR,"Illegal special_bonds command");
      special_lj[1] = 0.0;
      special_lj[2] = 1.0;
      special_lj[3] = 1.0;
      special_coul[1] = 0.0;
      special_coul[2] = 1.0;
      special_coul[3] = 1.0;
      iarg += 1;
    } else if (strcmp(arg[iarg],"lj/coul") == 0) {
      if (iarg+4 > narg) error->all(FLERR,"Illegal special_bonds command");
      special_lj[1] = special_coul[1] = atof(arg[iarg+1]);
      special_lj[2] = special_coul[2] = atof(arg[iarg+2]);
      special_lj[3] = special_coul[3] = atof(arg[iarg+3]);
      iarg += 4;
    } else if (strcmp(arg[iarg],"lj") == 0) {
      if (iarg+4 > narg) error->all(FLERR,"Illegal special_bonds command");
      special_lj[1] = atof(arg[iarg+1]);
      special_lj[2] = atof(arg[iarg+2]);
      special_lj[3] = atof(arg[iarg+3]);
      iarg += 4;
    } else if (strcmp(arg[iarg],"coul") == 0) {
      if (iarg+4 > narg) error->all(FLERR,"Illegal special_bonds command");
      special_coul[1] = atof(arg[iarg+1]);
      special_coul[2] = atof(arg[iarg+2]);
      special_coul[3] = atof(arg[iarg+3]);
      iarg += 4;
    } else if (strcmp(arg[iarg],"angle") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal special_bonds command");
      if (strcmp(arg[iarg+1],"no") == 0) special_angle = 0;
      else if (strcmp(arg[iarg+1],"yes") == 0) special_angle = 1;
      else error->all(FLERR,"Illegal special_bonds command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"dihedral") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal special_bonds command");
      if (strcmp(arg[iarg+1],"no") == 0) special_dihedral = 0;
      else if (strcmp(arg[iarg+1],"yes") == 0) special_dihedral = 1;
      else error->all(FLERR,"Illegal special_bonds command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"extra") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal special_bonds command");
      special_extra = atoi(arg[iarg+1]);
      iarg += 2;
    } else error->all(FLERR,"Illegal special_bonds command");
  }

  for (int i = 1; i <= 3; i++)
    if (special_lj[i] < 0.0 || special_lj[i] > 1.0 ||
        special_coul[i] < 0.0 || special_coul[i] > 1.0)
      error->all(FLERR,"Illegal special_bonds command");

  if (special_extra < 0) error->all(FLERR,"Illegal special_bonds command");
}

/* ----------------------------------------------------------------------
   compute bounds implied by numeric str with a possible wildcard asterik
   1 = lower bound, nmax = upper bound
   5 possibilities:
     (1) i = i to i, (2) * = nmin to nmax,
     (3) i* = i to nmax, (4) *j = nmin to j, (5) i*j = i to j
   return nlo,nhi
------------------------------------------------------------------------- */

void Force::bounds(char *str, int nmax, int &nlo, int &nhi, int nmin)
{
  char *ptr = strchr(str,'*');

  if (ptr == NULL) {
    nlo = nhi = atoi(str);
  } else if (strlen(str) == 1) {
    nlo = nmin;
    nhi = nmax;
  } else if (ptr == str) {
    nlo = nmin;
    nhi = atoi(ptr+1);
  } else if (strlen(ptr+1) == 0) {
    nlo = atoi(str);
    nhi = nmax;
  } else {
    nlo = atoi(str);
    nhi = atoi(ptr+1);
  }

  if (nlo < nmin || nhi > nmax)
    error->all(FLERR,"Numeric index is out of bounds");
}

/* ----------------------------------------------------------------------
   read a floating point value from a string
   generate an error if not a legitimate floating point value
   called by various commands to check validity of their arguments
------------------------------------------------------------------------- */

double Force::numeric(const char *file, const int line, const char *const str)
{
    const unsigned int n = strlen(str);
    char *dstr;
    dstr = new char[n+1];
    for (unsigned int i = 0; i < n; i++)
    {
        if (isdigit(str[i]) ||
            str[i] == '-'   ||
            str[i] == '+'   ||
            str[i] == '.'   ||
            str[i] == 'e'   ||
            str[i] == 'E')
            dstr[i] = str[i];
        else if (str[i] == '\r' && i == n-1)
            dstr[i] = '\0';
        else
            error->all(file, line, "Expected floating point parameter in input script or data file");
    }
    dstr[n] = '\0';

    const double val = atof(dstr);
    delete [] dstr;
    return val;
}

/* ----------------------------------------------------------------------
   read an integer value from a string
   generate an error if not a legitimate integer value
   called by various commands to check validity of their arguments
------------------------------------------------------------------------- */

int Force::inumeric(const char *file, const int line, const char *const str)
{
    const unsigned int n = strlen(str);
    char *istr;
    istr = new char[n+1];
    for (unsigned int i = 0; i < n; i++)
    {
        if (isdigit(str[i]) || str[i] == '-' || str[i] == '+')
            istr[i] = str[i];
        else if (str[i] == '\r' && i == n-1)
            istr[i] = '\0';
        else
            error->all(file, line, "Expected integer parameter in input script or data file");
    }
    istr[n] = '\0';

    const int val = atoi(istr);
    delete [] istr;
    return val;
}

/* ----------------------------------------------------------------------
   memory usage of force classes
------------------------------------------------------------------------- */

bigint Force::memory_usage()
{
  bigint bytes = 0;
  if (pair) bytes += static_cast<bigint> (pair->memory_usage());
  if (bond) bytes += static_cast<bigint> (bond->memory_usage());
  if (angle) bytes += static_cast<bigint> (angle->memory_usage());
  if (dihedral) bytes += static_cast<bigint> (dihedral->memory_usage());
  if (improper) bytes += static_cast<bigint> (improper->memory_usage());
  if (kspace) bytes += static_cast<bigint> (kspace->memory_usage());
  return bytes;
}
