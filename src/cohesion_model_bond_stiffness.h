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
    Valentin Baric (University of Bremen, Germany)
    v.baric@iwt.uni-bremen.de

    Copyright 2014-     DCS Computing GmbH, Linz
------------------------------------------------------------------------- */

#ifdef COHESION_MODEL
COHESION_MODEL(COHESION_BOND_STIFFNESS,bond/stiffness,10)
#else

#ifndef COHESION_MODEL_BOND_STIFFNESS_H_
#define COHESION_MODEL_BOND_STIFFNESS_H_

#include "contact_models.h"
#include "math.h"
#include "math_extra_liggghts.h"
#include "global_properties.h"
#include "fix_property_atom.h"
#include "neighbor.h"
#include "pointers.h"
#include <stdio.h>
#include <string>
#define MAXLINE 256

namespace MODEL_PARAMS
{
    inline static MatrixProperty* createBondTensileStrength2(PropertyRegistry & registry, const char * caller, bool sanity_checks)
    {
      return createPerTypePairProperty(registry, "bondTensileStrength", caller);
    }

    inline static MatrixProperty* createBondShearStrength2(PropertyRegistry & registry, const char * caller, bool sanity_checks)
    {
      return createPerTypePairProperty(registry, "bondShearStrength", caller);
    }

    inline static MatrixProperty* createBondRadius2(PropertyRegistry & registry, const char * caller, bool sanity_checks)
    {
      return createPerTypePairProperty(registry, "bondRadius", caller);
    }
    inline static MatrixProperty* createBondMaxSeparationDistance2(PropertyRegistry & registry, const char * caller, bool sanity_checks)
    {
     return createPerTypePairProperty(registry, "bondMaxSeparationDistance", caller);
    }
    inline static MatrixProperty* createBondSn(PropertyRegistry & registry, const char * caller, bool sanity_checks)
    {
     return createPerTypePairProperty(registry, "bondSn", caller);
    }
    inline static MatrixProperty* createBondSt(PropertyRegistry & registry, const char * caller, bool sanity_checks)
    {
     return createPerTypePairProperty(registry, "bondSt", caller);
    }
    inline static MatrixProperty* createBondRf(PropertyRegistry & registry, const char * caller, bool sanity_checks)
    {
     return createPerTypePairProperty(registry, "bondRf", caller);
    }
}

namespace LIGGGHTS {

namespace ContactModels {

  template<>
  class CohesionModel<COHESION_BOND_STIFFNESS> : public CohesionModelBase {

  public:
    CohesionModel(LAMMPS * lmp, IContactHistorySetup * hsetup, class ContactModelBase * cmb) :
        CohesionModelBase(lmp, hsetup, cmb),
        bondTensileStrength(NULL), 
        bondShearStrength(NULL),
        bondMaxSeparationDistance(NULL), 
        history_contact(0), 
        bondSn(NULL), 
        bondSt(NULL), 
        bondRf(NULL)
    {
      history_contact = hsetup->add_history_value("bondhistnormforcex","0");
      hsetup->add_history_value("bondhistnormforcey","0");
      hsetup->add_history_value("bondhistnormforcez","0");
      hsetup->add_history_value("bondhisttanforcex","0");
      hsetup->add_history_value("bondhisttanforcey","0");
      hsetup->add_history_value("bondhisttanforcez","0");
      hsetup->add_history_value("bondhistnormtorquex","0");
      hsetup->add_history_value("bondhistnormtorquey","0");
      hsetup->add_history_value("bondhistnormtorquez","0");
      hsetup->add_history_value("bondhisttantorquex","0");
      hsetup->add_history_value("bondhisttantorquey","0");
      hsetup->add_history_value("bondhisttantorquez","0");
      hsetup->add_history_value("contflag","0");
      hsetup->add_history_value("switchflag","0");
       if(cmb->is_wall())
        error->warning(FLERR,"The cohesion model bond/stiffness does not support walls");
    }

    void registerSettings(Settings& settings)
    {
      settings.registerOnOff("per_molecule",perMolecule_,true);
      settings.registerOnOff("full_length",fullLength_,true);
      settings.registerOnOff("write_file",writeFile_,true);
    }

    inline void postSettings(IContactHistorySetup * hsetup, ContactModelBase *cmb) {}

    void connectToProperties(PropertyRegistry & registry) {
      registry.registerProperty("bondTensileStrength", &MODEL_PARAMS::createBondTensileStrength2);
      registry.registerProperty("bondShearStrength", &MODEL_PARAMS::createBondShearStrength2);
      registry.registerProperty("bondMaxSeparationDistance", &MODEL_PARAMS::createBondMaxSeparationDistance2);
      registry.registerProperty("bondSn", &MODEL_PARAMS::createBondSn);
      registry.registerProperty("bondSt", &MODEL_PARAMS::createBondSt);
      registry.registerProperty("bondRf", &MODEL_PARAMS::createBondRf);

      registry.connect("bondTensileStrength", bondTensileStrength,"cohesion_model bond/stiffness");
      registry.connect("bondShearStrength", bondShearStrength,"cohesion_model bond/stiffness");
      registry.connect("bondMaxSeparationDistance", bondMaxSeparationDistance,"cohesion_model bond/stiffness");
      registry.connect("bondSn", bondSn,"cohesion_model bond/stiffness");
      registry.connect("bondSt", bondSt,"cohesion_model bond/stiffness");
      registry.connect("bondRf", bondRf,"cohesion_model bond/stiffness");

      // error checks on coarsegraining
      if(force->cg_active())
        error->cg(FLERR,"cohesion model bond/stiffness");

      const char* neigharg[2];
      neigharg[0] = "contact_distance_factor";
      char arg2[30];
      sprintf(arg2,"%e", 1.2);
      neigharg[1] = arg2;
      neighbor->modify_params(2,const_cast<char**>(neigharg));
    }
    inline void endSurfacesIntersect(SurfacesIntersectData &sidata, ForceData&, ForceData&) {}
    void beginPass(SurfacesIntersectData&, ForceData&, ForceData&){}
    void endPass(SurfacesIntersectData&, ForceData&, ForceData&){}

    float calc_phi(double dx, double dy, double dz)
    {
    // Calculate the angle phi in spherical coordinates
    // Depending on the distance between the point dx, dy, dz
      if(dx > 0) return atan(dx/dy);
      else if(dx == 0)
      {
        if (dy > 0) return M_PI / 2;
        else return - M_PI / 2;
      }
      else if(dx < 0 && dy >= 0) return atan(dy/dx) + M_PI;
      else return atan(dy/dx) - M_PI;
    }


    void surfacesIntersect(SurfacesIntersectData & sidata, ForceData & i_forces, ForceData & j_forces)
    {
      /*
      If the surfaces of the particles are intersecting there are three options that can be applied:

      1.  We are in timestep < 1 and bonds can be set. This accounts for slight fluctuations
          in distance during the initialization. This option will not be called
          after this time. If the option per_molecule is selected, bonds are only set
          between particles of the same molecule
      2.  Timestep > 1 and there is no bond. Nothing happens, no interaction by bonds
      3.  Timestep > 1 and an active bond. The force balances will be applied
          and if the breaking criteria is fulfilled, the bond breaks
      */
      // If the bonds are only to be set within molecules the rest is skipped here
      if (perMolecule_ && !(atom->molecule[sidata.i] == atom->molecule[sidata.j])){
        if(sidata.contact_flags) *sidata.contact_flags &= ~CONTACT_COHESION_MODEL;
        return;
      }

      int i = sidata.i;
      int j = sidata.j;
      double radi = sidata.radi;
      double radj = sidata.radj;
      // Particle coordinates
      double **x = atom->x;
      // Particle velocities
      double **v = atom->v;
      double **omega = atom->omega;

           // Timestep
      double dt = update->dt;

      // Vectors for normal and tangential force
      double dnforce[3],dtforce[3];
      // Vectors for normal and tangential torque
      double dntorque[3],dttorque[3];
      // Rotation of the coordinate system
      double rot;
      // boolean for an active bond and the destruction of a bond
      bool bond_active = false;

      if(sidata.contact_flags) *sidata.contact_flags |= CONTACT_COHESION_MODEL;
      // The history values for the forces and contactflag
      double * const historyvalues = &sidata.contact_history[history_contact];
      // The values for the bond radius
      double rbond = bondRf[sidata.itype][sidata.jtype] * MIN(radi,radj);
      // Values for the bond stiffness in normal and tangential direction
      double Sn = bondSn[sidata.itype][sidata.jtype];
      double St = bondSt[sidata.itype][sidata.jtype];
      // Bond values for tensile strength
      double sigman = bondTensileStrength[sidata.itype][sidata.jtype];
      double sigmat = bondShearStrength[sidata.itype][sidata.jtype];
      // Maximum distance of the particles that can be connected by a bond
      double maxDistance = bondMaxSeparationDistance[sidata.itype][sidata.jtype];

      bool switchid = false;

      // Calculation of bond properties, such as length and stretching
      // The area section of the bond
      //////////////////////////////////////////////////////////////////
      // 2. criteria: No active bond and timestep > 100:
      //              Cancel here
      if((update->ntimestep > 1) && !MathExtraLiggghts::compDouble(historyvalues[12],1.0,1e-6))
      {
        if(sidata.contact_flags) *sidata.contact_flags &= ~CONTACT_COHESION_MODEL;
        return;
      }
      if ((update->ntimestep > 0) & MathExtraLiggghts::compDouble(historyvalues[12],1.0,1e-6))
      {
        if ((historyvalues[13] == 0 && atom->tag[sidata.i] > atom->tag[sidata.j]) ||
            (historyvalues[13] == 1 && atom->tag[sidata.i] < atom->tag[sidata.j]))
        {
          switchid = true;
          // Particle velocities
          // Index of the first particle
          i = sidata.j;
          // Index of the second particle
          j = sidata.i;
          // Radius of the first particle
          radi = sidata.radj;
          // Radius of the second particle
          radj = sidata.radi;
        }else{
           // Particle velocities
           // Index of the first particle
          i = sidata.i;
          // Index of the second particle
          j = sidata.j;
          // Radius of the first particle
          radi = sidata.radi;
          // Radius of the second particle
          radj = sidata.radj;
        }
      }
      //////////////////////////////////////////////////////////////////

      const double A = M_PI * rbond * rbond;
      // The Inertia Value
      const double J = A * 0.5 * rbond * rbond;
      double dx = x[i][0] - x[j][0];
      double dy = x[i][1] - x[j][1];
      double dz = x[i][2] - x[j][2];
      const double dist = sqrt(dx*dx + dy*dy + dz*dz);

            // If fullLength_ is activated the distance of the two centers of mass
      // of the particles is used for the bond length
      // If it is deactivated the bond starts where the skin of the bond touches
      // the respective particle surface

      if(!fullLength_)
      {
        // The deviation of the bond length from the
        // distance from the center of mass
        /*
        double k1 = sqrt(radi*radi - rbond*rbond);
        double k2 = sqrt(radj*radj - rbond*rbond);
        // The angles between the two particles
        double phi1 = calc_phi(dx,dy,dz);
        double phi2 = calc_phi(-dx,-dy,-dz);
        double theta1 = acos(dz / dist);
        double theta2 = acos(-dz / dist);
        // Rotated length of the cut sections
        double dxr1 = k1 * sin(theta1) * cos(phi1);
        double dyr1 = k1 * sin(theta1) * sin(phi1);
        double dzr1 = k1 * cos(theta1);
        double dxr2 = k2 * sin(theta2) * cos(phi2);
        double dyr2 = k2 * sin(theta2) * sin(phi2);
        double dzr2 = k2 * cos(theta2);
        */
        double bLen = dist - radi - radj;
        // Real bond length
        dx = bLen * dx / dist;
        dy = bLen * dy / dist;
        dz = bLen * dz / dist;
      }
      // The length of the bond
      const double r = sqrt(dx * dx + dy * dy + dz * dz);
      // Inverted length
      const double rinv = 1./r;
      // squared length
      const double rsq = r * r;
      // inverted squared length
      const double rsqinv = 1./rsq;

      /////////////////////////////////////////////////////////////////
      // 1. criteria: Definition of a bond
      if ((dist < maxDistance) && (update->ntimestep < 1))
      {
        // The initial index distribution has to be stored
        // Because if the atom list is reordered later, it might switch i and j.
        // This would lead to the situation, that forces are initially added to the first atom
        // but later be subtracted or vice versa
        if (atom->tag[i] < atom->tag[j])
        {
          historyvalues[13] = 0;
        }
        else
        {
          historyvalues[13] = 1;
        }

        historyvalues[12] = 1.0;
        *sidata.contact_flags |= CONTACT_COHESION_MODEL;
        bond_active = true;
        fprintf(screen, "Generating Bond at timestep " BIGINT_FORMAT " between %d and %d\n", update->ntimestep, atom->tag[i], atom->tag[j]);
        if (writeFile_){
          FILE * pfile;
          pfile = fopen ("generated_bonds.txt","a");
          if (pfile != NULL)
            {
              fprintf(pfile,BIGINT_FORMAT " %d %d %e %e %e\n", update->ntimestep, atom->tag[i], atom->tag[j],radi,radj,r);
              fclose (pfile);
            }
        }
      }
      // 3. criteria, existing bond: calculate it
      else if (MathExtraLiggghts::compDouble(historyvalues[12],1.0,1e-6))
      {
        bond_active = true;
      }

      //////////////////////////////////////////////////////////////////
      //////////////////////////////////////////////////////////////////

      if(bond_active)
      {
      //relative translational velocity
      const double vr1 = v[i][0] - v[j][0];
      const double vr2 = v[i][1] - v[j][1];
      const double vr3 = v[i][2] - v[j][2];


      const double enx = dx * rinv;
      const double eny = dy * rinv;
      const double enz = dz * rinv;

      // normal component
      const double vn = vr1 * enx + vr2 * eny + vr3 * enz;
      const double vn1 = vn * enx;
      const double vn2 = vn * eny;
      const double vn3 = vn * enz;

      // tangential component
      const double vt1 = vr1 - vn1;
      const double vt2 = vr2 - vn2;
      const double vt3 = vr3 - vn3;

      // relative rotational velocity for shear
      double wr1 = (radi * omega[i][0] + radj * omega[j][0]) * rinv;
      double wr2 = (radi * omega[i][1] + radj * omega[j][1]) * rinv;
      double wr3 = (radi * omega[i][2] + radj * omega[j][2]) * rinv;

      // relative velocities for shear
      const double vtr1 = vt1 - (dz * wr2 - dy * wr3);
      const double vtr2 = vt2 - (dx * wr3 - dz * wr1);
      const double vtr3 = vt3 - (dy * wr1 - dx * wr2);

      // relative rotational velocity for torsion and bending
      wr1 = (radi * omega[i][0] - radj * omega[j][0]) * rinv;
      wr2 = (radi * omega[i][1] - radj * omega[j][1]) * rinv;
      wr3 = (radi * omega[i][2] - radj * omega[j][2]) * rinv;

      // normal component
      const double wnnr = wr1 * dx + wr2 * dy + wr3 * dz;
      const double wn1 = dx * wnnr * rsqinv;
      const double wn2 = dy * wnnr * rsqinv;
      const double wn3 = dz * wnnr * rsqinv;

      // tangential component
      const double wt1 = wr1 - wn1;
      const double wt2 = wr2 - wn2;
      const double wt3 = wr3 - wn3;

      // calc change in normal forces
      dnforce[0] = - vn1 * Sn * A * dt;
      dnforce[1] = - vn2 * Sn * A * dt;
      dnforce[2] = - vn3 * Sn * A * dt;

      // calc change in shear forces
      dtforce[0] = - vtr1 * St * A * dt;
      dtforce[1] = - vtr2 * St * A * dt;
      dtforce[2] = - vtr3 * St * A * dt;

      // calc change in normal torque
      dntorque[0] = - wn1 * St * J * dt;
      dntorque[1] = - wn2 * St * J * dt;
      dntorque[2] = - wn3 * St * J * dt;

      // calc change in tang torque
      dttorque[0] = - wt1 * Sn * J * 0.5 * dt;
      dttorque[1] = - wt2 * Sn * J * 0.5 * dt;
      dttorque[2] = - wt3 * Sn * J * 0.5 * dt;

      // rotate forces
      //rotate normal force
      rot = historyvalues[0] * dx + historyvalues[1] * dy + historyvalues[2] * dz;
      rot *= rsqinv;
      historyvalues[0] = rot * dx;
      historyvalues[1] = rot * dy;
      historyvalues[2] = rot * dz;

      //rotate tangential force
      rot = historyvalues[3]*dx + historyvalues[4]*dy + historyvalues[5]*dz;
      rot *= rsqinv;
      historyvalues[3] -= rot * dx;
      historyvalues[4] -= rot * dy;
      historyvalues[5] -= rot * dz;

      //rotate normal torque
      rot = historyvalues[6]*dx + historyvalues[7]*dy + historyvalues[8]*dz;
      rot *= rsqinv;
      historyvalues[6] = rot * dx;
      historyvalues[7] = rot * dy;
      historyvalues[8] = rot * dz;

      //rotate tangential torque
      rot = historyvalues[9]*dx + historyvalues[10]*dy + historyvalues[11]*dz;
      rot *= rsqinv;
      historyvalues[ 9] -= rot * dx;
      historyvalues[10] -= rot * dy;
      historyvalues[11] -= rot * dz;

      //increment normal and tangential force and torque
      historyvalues[ 0] += dnforce[0];
      historyvalues[ 1] += dnforce[1];
      historyvalues[ 2] += dnforce[2];
      historyvalues[ 3] += dtforce[0];
      historyvalues[ 4] += dtforce[1];
      historyvalues[ 5] += dtforce[2];
      historyvalues[ 6] += dntorque[0];
      historyvalues[ 7] += dntorque[1];
      historyvalues[ 8] += dntorque[2];
      historyvalues[ 9] += dttorque[0];
      historyvalues[10] += dttorque[1];
      historyvalues[11] += dttorque[2];

      const double tor1 = - rinv * (dy * historyvalues[5] - dz * historyvalues[4]);
      const double tor2 = - rinv * (dz * historyvalues[3] - dx * historyvalues[5]);
      const double tor3 = - rinv * (dx * historyvalues[4] - dy * historyvalues[3]);

      // The effective normal, tangential force and torque
      double nforce_mag = vectorMag3D(&historyvalues[0]);
      double tforce_mag = vectorMag3D(&historyvalues[3]);
      double ntorque_mag = vectorMag3D(&historyvalues[6]);
      double ttorque_mag = vectorMag3D(&historyvalues[9]);

      // Boolean values if a bond is broken
      bool nstress = sigman < (nforce_mag/A + 2.*ttorque_mag/J*rbond);
      bool tstress = sigmat    < (tforce_mag/A +    ntorque_mag/J*rbond);

      if(nstress || tstress)
      {
        if(sidata.contact_flags) *sidata.contact_flags &= ~CONTACT_COHESION_MODEL;
        historyvalues[12] = 0;
        fprintf(screen, "broken bond between atoms %d and %d at step " BIGINT_FORMAT,atom->tag[i],atom->tag[j],update->ntimestep);
          if(nstress)fprintf(screen," (nstress) %e", (nforce_mag/A + 2.*ttorque_mag/J*rbond));
          if(tstress)fprintf(screen," (tstress) %e",(tforce_mag/A +    ntorque_mag/J*rbond));
        fprintf(screen,"\n");

        if (writeFile_){
          FILE * pfile;
          pfile = fopen("broken_bonds.txt","a");

          if (pfile != NULL)
          {
            if (nstress && !tstress){
              fprintf (pfile, BIGINT_FORMAT " %d %d nstress\n", update->ntimestep,atom->tag[i],atom->tag[j]);
            }
            else if (!nstress && tstress){
              fprintf (pfile, BIGINT_FORMAT " %d %d tstress\n", update->ntimestep,atom->tag[i],atom->tag[j]);
            }
            else
            {
              fprintf (pfile, BIGINT_FORMAT " %d %d tstress + nstress\n", update->ntimestep,atom->tag[i],atom->tag[j]);
            }
            fclose(pfile);
          }
        }
      }
      sidata.has_force_update = true;
      // return resulting forces
      if(sidata.is_wall) {
        // Not implemented
      } else {
          if (!switchid)
          {
            i_forces.delta_F[0] += (historyvalues[0] + historyvalues[3]);
            i_forces.delta_F[1] += (historyvalues[1] + historyvalues[4]);
            i_forces.delta_F[2] += (historyvalues[2] + historyvalues[5]);
            i_forces.delta_torque[0] += radi * tor1 + (historyvalues[6] + historyvalues[ 9]);
            i_forces.delta_torque[1] += radi * tor2 + (historyvalues[7] + historyvalues[10]);
            i_forces.delta_torque[2] += radi * tor3 + (historyvalues[8] + historyvalues[11]);

            j_forces.delta_F[0] -= (historyvalues[0] + historyvalues[3]);
            j_forces.delta_F[1] -= (historyvalues[1] + historyvalues[4]);
            j_forces.delta_F[2] -= (historyvalues[2] + historyvalues[5]);
            j_forces.delta_torque[0] += radj * tor1 - (historyvalues[6] + historyvalues[ 9]);
            j_forces.delta_torque[1] += radj * tor2 - (historyvalues[7] + historyvalues[10]);
            j_forces.delta_torque[2] += radj * tor3 - (historyvalues[8] + historyvalues[11]);
          }else{
            j_forces.delta_F[0] += (historyvalues[0] + historyvalues[3]);
            j_forces.delta_F[1] += (historyvalues[1] + historyvalues[4]);
            j_forces.delta_F[2] += (historyvalues[2] + historyvalues[5]);
            j_forces.delta_torque[0] += radi * tor1 + (historyvalues[6] + historyvalues[ 9]);
            j_forces.delta_torque[1] += radi * tor2 + (historyvalues[7] + historyvalues[10]);
            j_forces.delta_torque[2] += radi * tor3 + (historyvalues[8] + historyvalues[11]);

            i_forces.delta_F[0] -= (historyvalues[0] + historyvalues[3]);
            i_forces.delta_F[1] -= (historyvalues[1] + historyvalues[4]);
            i_forces.delta_F[2] -= (historyvalues[2] + historyvalues[5]);
            i_forces.delta_torque[0] += radj * tor1 - (historyvalues[6] + historyvalues[ 9]);
            i_forces.delta_torque[1] += radj * tor2 - (historyvalues[7] + historyvalues[10]);
            i_forces.delta_torque[2] += radj * tor3 - (historyvalues[8] + historyvalues[11]);
          }
        }
    }
    }

    void surfacesClose(SurfacesCloseData & scdata, ForceData & i_forces, ForceData & j_forces)
    {
      /*
      If the surfaces of the aprticles are close to each other there are three options that can be applied:

      1.  We are in timestep < 10 and bonds can be set. This accounts for slight fluctuations
          in distance during the initialization. This option will not be called
          after this time. If the option per_molecule is selected, bonds are only set
          between particles of the same molecule
      2.  Timestep > 1 and there is no bond. Nothing happens, no interaction by bonds
      3.  Timestep > 1 and an active bond. The force balances will be applied
          and if the breaking criteria is fulfilled, the bond breaks
      */
      // If the bonds are only to be set within molecules the rest is skipped here
      if (perMolecule_ && !(atom->molecule[scdata.i] == atom->molecule[scdata.j])){
        if(scdata.contact_flags) *scdata.contact_flags &= ~CONTACT_COHESION_MODEL;
        return;
      }
      int i = scdata.i;
      int j = scdata.j;
      double radi = scdata.radi;
      double radj = scdata.radj;
      // Timestep
      double dt = update->dt;
      // Vectors for normal and tangential force
      double dnforce[3],dtforce[3];
      // Vectors for normal and tangential torque
      double dntorque[3],dttorque[3];
      // Rotation of the coordinate system
      double rot;
      // boolean for an active bond and the destruction of a bond
      bool bond_active = false;

      if(scdata.contact_flags) *scdata.contact_flags |= CONTACT_COHESION_MODEL;
      // The history value for an active contact

      double * const historyvalues = &scdata.contact_history[history_contact];
      // The values for the bond radius
      double rbond = bondRf[scdata.itype][scdata.jtype] * MIN(radi,radj);
      // Bond stiffness in nromal and tangential direction
      double Sn = bondSn[scdata.itype][scdata.jtype];
      double St = bondSt[scdata.itype][scdata.jtype];
      // Bond values for tensile strength
      double sigman = bondTensileStrength[scdata.itype][scdata.jtype];
      double sigmat = bondShearStrength[scdata.itype][scdata.jtype];
      // Maximum distance of the particles that can be connected by a bond
      double maxDistance = bondMaxSeparationDistance[scdata.itype][scdata.jtype];

      bool switchid = false;
      // Particle coordinates
      double **x = atom->x;
      // particle velocities
      double **v = atom->v;
      double **omega = atom->omega;

      // Calculation of bond properties, such as length and stretching
      // The area section of the bond
      //////////////////////////////////////////////////////////////////
      // 2. criteria: No active bond and timestep > 100:
      //              Cancel here
      if((update->ntimestep > 0) && !MathExtraLiggghts::compDouble(historyvalues[12],1.0,1e-6))
      {
        if(scdata.contact_flags) *scdata.contact_flags &= ~CONTACT_COHESION_MODEL;
        return;
      }
      if ((update->ntimestep > 0) & MathExtraLiggghts::compDouble(historyvalues[12],1.0,1e-6))
      {
        if ((historyvalues[13] == 0 && atom->tag[scdata.i] > atom->tag[scdata.j]) ||
            (historyvalues[13] == 1 && atom->tag[scdata.i] < atom->tag[scdata.j]))
        {
          switchid = true;
          // Particle velocities
          // Index of the first particle
          i = scdata.j;
          // Index of the second particle
          j = scdata.i;
          // Radius of the first particle
          radi = scdata.radj;
          // Radius of the second particle
          radj = scdata.radi;
        }else{
           // Particle velocities
           // Index of the first particle
          i = scdata.i;
          // Index of the second particle
          j = scdata.j;
          // Radius of the first particle
          radi = scdata.radi;
          // Radius of the second particle
          radj = scdata.radj;
        }
      }
      //////////////////////////////////////////////////////////////////

      const double A = M_PI * rbond * rbond;
      // The Inertia Value
      const double J = A * 0.5 * rbond * rbond;
      double dx = x[i][0] - x[j][0];
      double dy = x[i][1] - x[j][1];
      double dz = x[i][2] - x[j][2];
      const double dist = sqrt(dx * dx + dy * dy + dz * dz);
      // If fullLength_ is activated the distance of the two centers of mass
      // of the particles is used for the bond length
      // If it is deactivated the bond starts where the skin of the bond touches
      // the respective particle surface
      if(!fullLength_)
      {
        // The deviation of the bond length from the
        // distance from the center of mass
        /*
        double k1 = sqrt(radi*radi - rbond*rbond);
        double k2 = sqrt(radj*radj - rbond*rbond);
        // The angles between the two particles
        double phi1 = calc_phi(dx,dy,dz);
        double phi2 = calc_phi(-dx,-dy,-dz);
        double theta1 = acos(dz / dist);
        double theta2 = acos(-dz / dist);
        // Rotated length of the cut sections
        double dxr1 = k1 * sin(theta1) * cos(phi1);
        double dyr1 = k1 * sin(theta1) * sin(phi1);
        double dzr1 = k1 * cos(theta1);
        double dxr2 = k2 * sin(theta2) * cos(phi2);
        double dyr2 = k2 * sin(theta2) * sin(phi2);
        double dzr2 = k2 * cos(theta2);
        */
        double bLen = dist - radi - radj;
        // Real bond length
        dx = bLen * dx / dist;
        dy = bLen * dy / dist;
        dz = bLen * dz / dist;
      }
      // The length of the bond
      // The length is count from the particle surfaces and not the center
      // of mass
      const double r = sqrt(dx * dx + dy * dy + dz * dz);
      // Inverted length
      const double rinv = 1./r;
      // squared length
      const double rsq = r*r;
      // inverted squared length
      const double rsqinv = 1./rsq;

      /////////////////////////////////////////////////////////////////
      // 1. criteria: Definition of a bond
      if ((dist < maxDistance) && (update->ntimestep == 0))
      {
        // The initial index distribution has to be stored
        // Because if the atom list is reordered later, it might switch i and j.
        // This would lead to the situation, that forces are initially added to the first atom
        // but later be subtracted or vice versa
        if (atom->tag[i] < atom->tag[j])
        {
          historyvalues[13] = 0;
        }
        else
        {
          historyvalues[13] = 1;
        }
      // If the bonds are only defined within an molecule
        historyvalues[12] = 1.0;
        *scdata.contact_flags |= CONTACT_COHESION_MODEL;
        bond_active = true;
        fprintf(screen, "Generating Bond at timestep " BIGINT_FORMAT " between %d and %d\n", update->ntimestep, atom->tag[i], atom->tag[j]);
        if (writeFile_){
          FILE * pfile;
          pfile = fopen ("generated_bonds.txt","a");
          if (pfile != NULL)
            {
              fprintf(pfile,BIGINT_FORMAT " %d %d %e %e %e\n", update->ntimestep, atom->tag[i], atom->tag[j],radi,radj,r);
              fclose (pfile);
            }
        }
      }
      // 3. criteria, existing bond: calculate it
      else if (MathExtraLiggghts::compDouble(historyvalues[12],1.0,1e-6))
      {
        bond_active = true;
      }
      //////////////////////////////////////////////////////////////////
      //////////////////////////////////////////////////////////////////
      // Disable calculation if there is no bond
      if(!MathExtraLiggghts::compDouble(historyvalues[12],1.0,1e-6) && (update->ntimestep > 1))
      {
          if(scdata.contact_flags) *scdata.contact_flags &= ~CONTACT_COHESION_MODEL;
          return;
      }
      //////////////////////////////////////////////////////////////////
      //////////////////////////////////////////////////////////////////
      // 3rd criteria: There is an existing bond so calculate its forces
      if(bond_active)
      {
      //relative translational velocity
      const double vr1 = v[i][0] - v[j][0];
      const double vr2 = v[i][1] - v[j][1];
      const double vr3 = v[i][2] - v[j][2];

      const double enx = dx * rinv;
      const double eny = dy * rinv;
      const double enz = dz * rinv;

      // normal component
      const double vn = vr1 * enx + vr2 * eny + vr3 * enz;
      const double vn1 = vn * enx;
      const double vn2 = vn * eny;
      const double vn3 = vn * enz;

      // tangential component
      const double vt1 = vr1 - vn1;
      const double vt2 = vr2 - vn2;
      const double vt3 = vr3 - vn3;

      // relative rotational velocity for shear
      double wr1 = (radi * omega[i][0] + radj * omega[j][0]) * rinv;
      double wr2 = (radi * omega[i][1] + radj * omega[j][1]) * rinv;
      double wr3 = (radi * omega[i][2] + radj * omega[j][2]) * rinv;

      // relative velocities for shear
      const double vtr1 = vt1 - (dz * wr2 - dy * wr3);
      const double vtr2 = vt2 - (dx * wr3 - dz * wr1);
      const double vtr3 = vt3 - (dy * wr1 - dx * wr2);

      // relative rotational velocity for torsion and bending
      wr1 = (radi * omega[i][0] - radj * omega[j][0]) * rinv;
      wr2 = (radi * omega[i][1] - radj * omega[j][1]) * rinv;
      wr3 = (radi * omega[i][2] - radj * omega[j][2]) * rinv;

      // normal component
      const double wn = wr1 * dx + wr2 * dy + wr3 * dz;
      const double wn1 = wn * enx;
      const double wn2 = wn * eny;
      const double wn3 = wn * enz;

      // tangential component
      const double wt1 = wr1 - wn1;
      const double wt2 = wr2 - wn2;
      const double wt3 = wr3 - wn3;

      // calc change in normal forces
      dnforce[0] = - vn1 * Sn * A * dt;
      dnforce[1] = - vn2 * Sn * A * dt;
      dnforce[2] = - vn3 * Sn * A * dt;

      // calc change in shear forces
      dtforce[0] = - vtr1 * St * A * dt;
      dtforce[1] = - vtr2 * St * A * dt;
      dtforce[2] = - vtr3 * St * A * dt;

      // calc change in normal torque
      dntorque[0] = - wn1 * St * J * dt;
      dntorque[1] = - wn2 * St * J * dt;
      dntorque[2] = - wn3 * St * J * dt;

      // calc change in tang torque
      dttorque[0] = - wt1 * Sn * J * 0.5 * dt;
      dttorque[1] = - wt2 * Sn * J * 0.5 * dt;
      dttorque[2] = - wt3 * Sn * J * 0.5 * dt;

      // rotate forces
      //rotate normal force
      rot = historyvalues[0]*dx + historyvalues[1]*dy + historyvalues[2]*dz;
      rot *= rsqinv;
      historyvalues[0] = rot * dx;
      historyvalues[1] = rot * dy;
      historyvalues[2] = rot * dz;

      //rotate tangential force
      rot = historyvalues[3]*dx + historyvalues[4]*dy + historyvalues[5]*dz;
      rot *= rsqinv;
      historyvalues[3] -= rot * dx;
      historyvalues[4] -= rot * dy;
      historyvalues[5] -= rot * dz;

      //rotate normal torque
      rot = historyvalues[6]*dx + historyvalues[7]*dy + historyvalues[8]*dz;
      rot *= rsqinv;
      historyvalues[6] = rot * dx;
      historyvalues[7] = rot * dy;
      historyvalues[8] = rot * dz;

      //rotate tangential torque
      rot = historyvalues[9]*dx + historyvalues[10]*dy + historyvalues[11]*dz;
      rot *= rsqinv;
      historyvalues[ 9] -= rot * dx;
      historyvalues[10] -= rot * dy;
      historyvalues[11] -= rot * dz;

      //increment normal and tangential force and torque
      historyvalues[ 0] += dnforce[0];
      historyvalues[ 1] += dnforce[1];
      historyvalues[ 2] += dnforce[2];
      historyvalues[ 3] += dtforce[0];
      historyvalues[ 4] += dtforce[1];
      historyvalues[ 5] += dtforce[2];
      historyvalues[ 6] += dntorque[0];
      historyvalues[ 7] += dntorque[1];
      historyvalues[ 8] += dntorque[2];
      historyvalues[ 9] += dttorque[0];
      historyvalues[10] += dttorque[1];
      historyvalues[11] += dttorque[2];

      const double tor1 = - rinv * (dy * historyvalues[5] - dz * historyvalues[4]);
      const double tor2 = - rinv * (dz * historyvalues[3] - dx * historyvalues[5]);
      const double tor3 = - rinv * (dx * historyvalues[4] - dy * historyvalues[3]);

      // The effective normal, tangential force and torque
      double nforce_mag = vectorMag3D(&historyvalues[0]);
      double tforce_mag = vectorMag3D(&historyvalues[3]);
      double ntorque_mag = vectorMag3D(&historyvalues[6]);
      double ttorque_mag = vectorMag3D(&historyvalues[9]);

      // Boolean values if a bond is broken
      bool nstress = sigman < (nforce_mag/A + 2.*ttorque_mag/J*rbond);
      bool tstress = sigmat    < (tforce_mag/A +    ntorque_mag/J*rbond);

      if(nstress || tstress)
      {
        historyvalues[12] = 0;
        fprintf(screen, "broken bond between atoms %d and %d at step " BIGINT_FORMAT,atom->tag[i],atom->tag[j],update->ntimestep);
          if(nstress)fprintf(screen," (nstress) %e", (nforce_mag/A + 2.*ttorque_mag/J*rbond));
          if(tstress)fprintf(screen," (tstress) %e",(tforce_mag/A +    ntorque_mag/J*rbond));
        fprintf(screen,"\n");
        if (writeFile_){
          FILE * pfile;
          pfile = fopen("broken_bonds.txt","a");
          if (pfile != NULL)
          {
            if (nstress && !tstress){
              fprintf (pfile,BIGINT_FORMAT " %d %d nstress\n", update->ntimestep,atom->tag[i],atom->tag[j]);
            }
            else if (!nstress && tstress){
              fprintf (pfile, BIGINT_FORMAT " %d %d tstress\n", update->ntimestep,atom->tag[i],atom->tag[j]);
            }
            else
            {
              fprintf (pfile, BIGINT_FORMAT " %d %d tstress + nstress\n", update->ntimestep,atom->tag[i],atom->tag[j]);
            }
            fclose(pfile);
          }
        }
        if(scdata.contact_flags) *scdata.contact_flags &= ~CONTACT_COHESION_MODEL;
        historyvalues[12] = 0.0;
      }
      scdata.has_force_update = true;
      // return resulting forces
      if(scdata.is_wall) {
        // Not implemented

      } else {
          if (!switchid)
          {
            i_forces.delta_F[0] += (historyvalues[0] + historyvalues[3]);
            i_forces.delta_F[1] += (historyvalues[1] + historyvalues[4]);
            i_forces.delta_F[2] += (historyvalues[2] + historyvalues[5]);
            i_forces.delta_torque[0] += radi * tor1 + (historyvalues[6] + historyvalues[ 9]);
            i_forces.delta_torque[1] += radi * tor2 + (historyvalues[7] + historyvalues[10]);
            i_forces.delta_torque[2] += radi * tor3 + (historyvalues[8] + historyvalues[11]);

            j_forces.delta_F[0] -= (historyvalues[0] + historyvalues[3]);
            j_forces.delta_F[1] -= (historyvalues[1] + historyvalues[4]);
            j_forces.delta_F[2] -= (historyvalues[2] + historyvalues[5]);
            j_forces.delta_torque[0] += radj * tor1 - (historyvalues[6] + historyvalues[ 9]);
            j_forces.delta_torque[1] += radj * tor2 - (historyvalues[7] + historyvalues[10]);
            j_forces.delta_torque[2] += radj * tor3 - (historyvalues[8] + historyvalues[11]);
          }else
          {
            j_forces.delta_F[0] += (historyvalues[0] + historyvalues[3]);
            j_forces.delta_F[1] += (historyvalues[1] + historyvalues[4]);
            j_forces.delta_F[2] += (historyvalues[2] + historyvalues[5]);
            j_forces.delta_torque[0] += radi * tor1 + (historyvalues[6] + historyvalues[ 9]);
            j_forces.delta_torque[1] += radi * tor2 + (historyvalues[7] + historyvalues[10]);
            j_forces.delta_torque[2] += radi * tor3 + (historyvalues[8] + historyvalues[11]);

            i_forces.delta_F[0] -= (historyvalues[0] + historyvalues[3]);
            i_forces.delta_F[1] -= (historyvalues[1] + historyvalues[4]);
            i_forces.delta_F[2] -= (historyvalues[2] + historyvalues[5]);
            i_forces.delta_torque[0] += radj * tor1 - (historyvalues[6] + historyvalues[ 9]);
            i_forces.delta_torque[1] += radj * tor2 - (historyvalues[7] + historyvalues[10]);
            i_forces.delta_torque[2] += radj * tor3 - (historyvalues[8] + historyvalues[11]);
          }

        }
    }
  }
  private:
    double ** bondTensileStrength;
    double ** bondShearStrength;
    double ** bondMaxSeparationDistance;
    double ** bondSn;
    double ** bondSt;
    double ** bondRf;
    int history_contact;
    bool perMolecule_,fullLength_, writeFile_;
  };
}
}
#endif
#endif
