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

    Copyright 2012-     DCS Computing GmbH, Linz
    Copyright 2009-2012 JKU Linz
------------------------------------------------------------------------- */

#include <string.h>
#include <stdlib.h>
#include "atom.h"
#include "update.h"
#include "respa.h"
#include "error.h"
#include "neighbor.h"
#include "memory.h"
#include "modify.h"
#include "group.h"
#include "comm.h"
#include <cmath>
#include "vector_liggghts.h"
#include "fix_cfd_coupling_convection.h"
#include "fix_property_atom.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixCfdCouplingConvection::FixCfdCouplingConvection(LAMMPS *lmp, int narg, char **arg) :
   Fix(lmp, narg, arg),
   is_convection_(true),
   fix_coupling_(0),
   fix_additionalFlux_(0),
   fix_heatFlux_(0)
{
    int iarg = 3;

    if(strstr(style,"radiation"))
        is_convection_ = false;

    if(narg < iarg + 2)
        error->fix_error(FLERR,this,"wrong number of arguments");
    if(strcmp(arg[iarg++],"T0") != 0)
        error->fix_error(FLERR,this,"expecting keyword 'T0'");
    T0_ = atof(arg[iarg++]);

    if(T0_ < 0.)
        error->fix_error(FLERR,this,"T0 must be >= 0");
}

/* ---------------------------------------------------------------------- */

FixCfdCouplingConvection::~FixCfdCouplingConvection()
{
}

/* ---------------------------------------------------------------------- */

void FixCfdCouplingConvection::pre_delete(bool unfixflag)
{
    if(fix_additionalFlux_ && is_convection_)
        modify->delete_fix("convectiveHeatFlux");
    if(fix_additionalFlux_ && !is_convection_)
        modify->delete_fix("radiativeHeatFlux");
}

/* ---------------------------------------------------------------------- */

int FixCfdCouplingConvection::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixCfdCouplingConvection::post_create()
{
  //  register convective/radiative flux
  if(!fix_additionalFlux_)
  {
        const char* fixarg[11];
        fixarg[0]=is_convection_?"convectiveHeatFlux":"radiativeHeatFlux";
        fixarg[1]="all";
        fixarg[2]="property/atom";
        fixarg[3]=is_convection_?"convectiveHeatFlux":"radiativeHeatFlux";
        fixarg[4]="scalar"; 
        fixarg[5]="no";    
        fixarg[6]="yes";    
        fixarg[7]="no";    
        fixarg[8]="0.";
        fix_additionalFlux_ = modify->add_fix_property_atom(9,const_cast<char**>(fixarg),style);
  }

  //  add heat transfer model if not yet active
  FixScalarTransportEquation *fix_ste = modify->find_fix_scalar_transport_equation("heattransfer");
  if(!fix_ste)
  {
        const char *newarg[15];
        newarg[0] = "ste_heattransfer";
        newarg[1] = group->names[igroup];
        newarg[2] = "transportequation/scalar";
        newarg[3] = "equation_id";
        newarg[4] = "heattransfer";
        newarg[5] = "quantity";
        newarg[6] = "Temp";
        newarg[7] = "default_value";
        char arg8[30];
        sprintf(arg8,"%f",T0_);
        newarg[8] = arg8;
        newarg[9] = "flux_quantity";
        newarg[10] = "heatFlux";
        newarg[11] = "source_quantity";
        newarg[12] = "heatSource";
        newarg[13] = "capacity_quantity";
        newarg[14] = "thermalCapacity";
        modify->add_fix(15,const_cast<char**>(newarg));
  }
}

/* ---------------------------------------------------------------------- */

void FixCfdCouplingConvection::init()
{
    // make sure there is only one fix of this style
    if(modify->n_fixes_style(style) != 1)
      error->fix_error(FLERR,this,"More than one fix of this style is not allowed");

    // find coupling fix
    fix_coupling_ = static_cast<FixCfdCoupling*>(modify->find_fix_style_strict("couple/cfd",0));
    if(!fix_coupling_)
      error->fix_error(FLERR,this,"needs a fix of type couple/cfd");

    //values to send to OF
    fix_coupling_->add_push_property("Temp","scalar-atom");

    //values to come from OF
    if(is_convection_)
        fix_coupling_->add_pull_property("convectiveHeatFlux","scalar-atom");
    else
        fix_coupling_->add_pull_property("radiativeHeatFlux","scalar-atom");

    // heat transfer added heatFlux, get reference to it
    fix_heatFlux_ = static_cast<FixPropertyAtom*>(modify->find_fix_property("heatFlux","property/atom","scalar",0,0,style));

    if(is_convection_)
        fix_additionalFlux_ = static_cast<FixPropertyAtom*>(modify->find_fix_property("convectiveHeatFlux","property/atom","scalar",0,0,style));
    else
        fix_additionalFlux_ = static_cast<FixPropertyAtom*>(modify->find_fix_property("radiativeHeatFlux","property/atom","scalar",0,0,style));
}

/* ---------------------------------------------------------------------- */

void FixCfdCouplingConvection::post_force(int)
{
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  // communicate convective/radiative flux to ghosts, there might be new data
  
  if(0 == neighbor->ago)
        fix_additionalFlux_->do_forward_comm();

  double *heatFlux = fix_heatFlux_->vector_atom;
  double *additionalFlux = fix_additionalFlux_->vector_atom;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit)
      heatFlux[i] += additionalFlux[i];
}
