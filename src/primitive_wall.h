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

    Christoph Kloss (DCS Computing GmbH, Linz)
    Christoph Kloss (JKU Linz)
    Philippe Seil (JKU Linz)

    Copyright 2012-     DCS Computing GmbH, Linz
    Copyright 2009-2012 JKU Linz
------------------------------------------------------------------------- */

#ifndef LMP_PRIMITIVE_WALL
#define LMP_PRIMITIVE_WALL

#include "container.h"
#include "neighbor.h"
#include "primitive_wall_definitions.h"

namespace LAMMPS_NS
{

  /*
   * class PrimitiveWall only holds logic for resolving contacts, neighbor list and contact tracking
   */

  class PrimitiveWall : protected Pointers
  {
      public:

        PrimitiveWall(LAMMPS *lmp,PRIMITIVE_WALL_DEFINITIONS::WallType wType_, int nParam_, double *param_)
        : Pointers(lmp), neighlist("neighlist"), wType(wType_), nParam(nParam_)
        {
            param = new double[nParam];
            for(int i=0;i<nParam;i++)
              param[i] = param_[i];
        }

        virtual ~PrimitiveWall()
        {
            delete []param;
        }

        inline int getNeighbors(int *&contactPtr);
        inline void handleContact(int iPart,double *&c_history, const double deltan);
        inline void handleNoContact(int iPart);
        inline void setContactHistorySize(int nPart);

        inline void buildNeighList(double neighCutoff, double **x, double *r, int nPart);

        inline double resolveContact(double *x, double r, double *delta);
        inline bool resolveNeighlist(double *x, double r, double treshold);

        inline int axis();
        inline double calcRadialDistance(double *pos, double *distvec);

        inline int isNear(int iPart,double treshold);

      private:
        ScalarContainer<int> neighlist;
        PRIMITIVE_WALL_DEFINITIONS::WallType wType;

        double *param;
        int nParam;

  };

  /*
   * implementation of class PrimitiveWall starts here
   */

  double PrimitiveWall::resolveContact(double *x, double r, double *delta)
  {
    return PRIMITIVE_WALL_DEFINITIONS::chooseContactTemplate(x, r, delta, param, wType);
  }

  int PrimitiveWall::axis()
  {
    return PRIMITIVE_WALL_DEFINITIONS::chooseAxis(wType);
  }

  double PrimitiveWall::calcRadialDistance(double *pos, double *distvec)
  {
    return PRIMITIVE_WALL_DEFINITIONS::chooseCalcRadialDistance(pos, param, distvec[0],distvec[1],distvec[2], wType);
  }

  bool PrimitiveWall::resolveNeighlist(double *x, double r, double treshold)
  {
    return PRIMITIVE_WALL_DEFINITIONS::chooseNeighlistTemplate(x,r,treshold,param,wType);
  }

  int PrimitiveWall::getNeighbors(int *&contactPtr)
  {
    contactPtr = neighlist.begin();
    return neighlist.size();
  }

  void PrimitiveWall::buildNeighList(double treshold, double **x, double *r, int nPart)
  {
    neighlist.clearContainer();
    for(int iPart=0;iPart<nPart;iPart++)
    {
      
      if(resolveNeighlist(x[iPart],r?r[iPart]:0.,treshold))
        neighlist.add(iPart);
    }
  }

  int PrimitiveWall::isNear(int iPart,double treshold)
  {
    if(resolveNeighlist(atom->x[iPart],atom->radius?atom->radius[iPart]:0.,treshold))
        return 1;

    return 0;
  }
} /* namespace LAMMPS_NS */
#endif /* PRIMITIVEWALL_H_ */
