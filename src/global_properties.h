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
    Richard Berger (JKU Linz)

    Copyright 2012-     DCS Computing GmbH, Linz
    Copyright 2009-2012 JKU Linz
------------------------------------------------------------------------- */

#ifndef GLOBAL_PROPERTIES_H_
#define GLOBAL_PROPERTIES_H_

#include "property_registry.h"

#define BIG 1.0e20

namespace MODEL_PARAMS
{

    /* -----------------------------------------------------------------------
     * Utility functions
     * ----------------------------------------------------------------------- */

    ScalarProperty* createScalarProperty(PropertyRegistry & registry, const char* name, const char * caller);
    ScalarProperty* createScalarProperty(PropertyRegistry & registry, const char* name, const char * caller, bool sanity_checks, const double lo = -BIG, const double hi = BIG);

    VectorProperty* createPerTypeProperty(PropertyRegistry & registry, const char* name, const char * caller);
    VectorProperty* createPerTypeProperty(PropertyRegistry & registry, const char* name, const char * caller, bool sanity_checks, const double lo = -BIG, const double hi = BIG);

    MatrixProperty* createPerTypePairProperty(PropertyRegistry & registry, const char * name, const char * caller);
    MatrixProperty* createPerTypePairProperty(PropertyRegistry & registry, const char * name, const char * caller, bool sanity_checks, const double lo = -BIG, const double hi = BIG);

    /* -----------------------------------------------------------------------
     * Property Creators
     * ----------------------------------------------------------------------- */
    VectorProperty* createCoalescenceModelSwitches(PropertyRegistry & registry, const char * caller, bool sanity_checks);
    VectorProperty* createCoalescenceModelSettings(PropertyRegistry & registry, const char * caller, bool sanity_checks);

    MatrixProperty* createCohesionEnergyDensity(PropertyRegistry & registry, const char * caller, bool sanity_checks);
    VectorProperty* createCohesionDistanceSettings(PropertyRegistry & registry, const char * caller, bool sanity_checks);
    VectorProperty* createCohesionModelSwitches(PropertyRegistry & registry, const char * caller, bool sanity_checks);
    ScalarProperty* createCharacteristicVelocity(PropertyRegistry & registry, const char * caller, bool sanity_checks);
    VectorProperty* createYoungsModulus(PropertyRegistry & registry, const char * caller, bool sanity_checks) ;
    VectorProperty* createPoissonsRatio(PropertyRegistry & registry, const char * caller, bool sanity_checks);
    MatrixProperty* createYeff(PropertyRegistry & registry, const char * caller, bool sanity_checks);
    MatrixProperty* createGeff(PropertyRegistry & registry, const char * caller, bool sanity_checks);
    MatrixProperty* createCoeffRest(PropertyRegistry & registry, const char * caller, bool sanity_checks);
    MatrixProperty* createCoeffRestLog(PropertyRegistry & registry, const char * caller, bool sanity_checks);
    MatrixProperty* createBetaEff(PropertyRegistry & registry, const char * caller, bool sanity_checks);

    MatrixProperty* createCoeffFrict(PropertyRegistry & registry, const char * caller, bool sanity_checks);
    MatrixProperty* createCoeffRollFrict(PropertyRegistry & registry, const char * caller, bool sanity_checks);
    MatrixProperty* createCoeffFricVisc(PropertyRegistry & registry, const char * caller, bool sanity_checks); // For new models
    MatrixProperty* createCoeffRollVisc(PropertyRegistry & registry, const char * caller, bool sanity_checks);
    MatrixProperty* createCoeffMu(PropertyRegistry & registry, const char * caller, bool sanity_checks);
    MatrixProperty* createCoeffRestMax(PropertyRegistry & registry, const char * caller, bool sanity_checks);
    MatrixProperty* createCoeffStc(PropertyRegistry & registry, const char * caller, bool sanity_checks);
    ScalarProperty* createRollingStiffness(PropertyRegistry & registry, const char * caller, bool sanity_checks);

    ScalarProperty* createLiquidVolume(PropertyRegistry & registry, const char * caller, bool sanity_checks);
    ScalarProperty* createLiquidDensity(PropertyRegistry & registry, const char * caller, bool sanity_checks);
    ScalarProperty* createSurfaceTension(PropertyRegistry & registry, const char * caller, bool sanity_checks);
    ScalarProperty* createSwitchModel(PropertyRegistry & registry, const char * caller, bool sanity_checks);
    ScalarProperty* createHistoryIndex(PropertyRegistry & registry, const char * caller, bool sanity_checks);
    VectorProperty* createContactAngle(PropertyRegistry & registry, const char * caller, bool sanity_checks);
    MatrixProperty* createKn(PropertyRegistry & registry, const char * caller, bool sanity_checks);
    MatrixProperty* createKt(PropertyRegistry & registry, const char * caller, bool sanity_checks);
    MatrixProperty* createGamman(PropertyRegistry & registry, const char * caller, bool sanity_checks);
    MatrixProperty* createGammat(PropertyRegistry & registry, const char * caller, bool sanity_checks);
    MatrixProperty* createGammanAbs(PropertyRegistry & registry, const char * caller, bool sanity_checks);
    MatrixProperty* createGammatAbs(PropertyRegistry & registry, const char * caller, bool sanity_checks);
    MatrixProperty* createCoeffMaxElasticStiffness(PropertyRegistry & registry, const char * caller, bool sanity_checks);
    MatrixProperty* createCoeffAdhesionStiffness(PropertyRegistry & registry, const char * caller, bool sanity_checks);
    MatrixProperty* createCoeffPlasticityDepth(PropertyRegistry & registry, const char * caller, bool sanity_checks);
    MatrixProperty* createPullOffForce(PropertyRegistry & registry, const char * caller, bool sanity_checks); // for luding's model
    MatrixProperty* createCoeffFrictionStiffness(PropertyRegistry & registry, const char * caller, bool sanity_checks); //for luding's tangential model
    MatrixProperty* createCoeffRollingStiffness(PropertyRegistry & registry, const char * caller, bool sanity_checks);
    ScalarProperty* createOverlapExponent(PropertyRegistry & registry, const char * caller, bool sanity_checks); //for Edinburgh's model
    ScalarProperty* createAdhesionExponent(PropertyRegistry & registry, const char * caller, bool sanity_checks);
    MatrixProperty* createSurfaceEnergy(PropertyRegistry & registry, const char * caller, bool sanity_checks);
    MatrixProperty* createTangentialMultiplier(PropertyRegistry & registry, const char * caller, bool sanity_checks);
    VectorProperty* createYieldRatio(PropertyRegistry & registry, const char * caller, bool sanity_checks); // thornton & ning
    MatrixProperty* createUnloadingStiffness(PropertyRegistry & registry, const char * caller, bool sanity_checks); // for hysteretic models 
    MatrixProperty* createLoadingStiffness(PropertyRegistry & registry, const char * caller, bool sanity_checks); // for hysteretic models 
    ScalarProperty* createRoughnessAbsolute(PropertyRegistry & registry, const char * caller, bool sanity_checks);
    ScalarProperty* createRoughnessRelative(PropertyRegistry & registry, const char * caller, bool sanity_checks);
    MatrixProperty* createNormalDampingCoefficient(PropertyRegistry & registry, const char * caller, bool sanity_checks);
    MatrixProperty* createTangentialDampingCoefficient(PropertyRegistry & registry, const char * caller, bool sanity_checks);
    MatrixProperty* createTangentialStiffness(PropertyRegistry & registry, const char * caller, bool sanity_checks);
    MatrixProperty* createInitialCohesiveStress(PropertyRegistry & registry, const char * caller, bool sanity_checks);
    MatrixProperty* createMaxCohesiveStress(PropertyRegistry & registry, const char * caller, bool sanity_checks);
    MatrixProperty* createCohesionStrength(PropertyRegistry & registry, const char * caller, bool sanity_checks);

    // fiber & bond model
    MatrixProperty* createDissipationMatrix(PropertyRegistry & registry, const char * name, const char * caller, bool sanity_checks);

    // Attrition models
    MatrixProperty* createAbrasiveWearSeverity(PropertyRegistry & registry, const char * caller, bool sanity_checks);
    MatrixProperty* createErosiveWearSeverity(PropertyRegistry & registry, const char * caller, bool sanity_checks);
    MatrixProperty* createVelocityLimit(PropertyRegistry & registry, const char * caller, bool sanity_checks);
    MatrixProperty* createForceLimit(PropertyRegistry & registry, const char * caller, bool sanity_checks);
    VectorProperty* createAttritionLimit(PropertyRegistry & registry, const char * caller, bool sanity_checks);
}

#endif /* GLOBAL_PROPERTIES_H_ */
