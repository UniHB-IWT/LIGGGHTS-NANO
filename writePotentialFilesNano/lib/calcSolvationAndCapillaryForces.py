# Import the library for calculating the capillary forces (separated for readability)
from lib.classCapillaryForce import capforce as calc_cap_force
# Numerical Python
import numpy as np
# System library
import sys
# Operating System library
import os
# Plotting
import matplotlib.pyplot as plt
# Integration for energy correspoinding to forces (required by LIGGGHTS)
from scipy import integrate
# Combine the particle radii to get all combinations
from itertools import combinations_with_replacement


class potential:
    """
        Calculate the interaction potentials between particles for DEM simulations.
    """

    def __init__(self, diameters,
                 tablefile='include.potential_table',
                 parameterfile='include.potential_params',
                 humidity=0.5,
                 equi_dist=0.57,
                 gamma_n=5e6,
                 gamma_t=0,
                 tau_s=0.65e9,
                 tau_r=0.25e9,
                 delta_h=0.15,
                 sigma_solv=0.21,
                 F_solv_b=-11.0,
                 F_solv_m=-2.4,
                 ):

        # Filename of the potentialfile (output)
        self.outputfile = tablefile
        # Filename of the potential parameterfile (output)
        self.parameterfile = parameterfile
        # Solvationforces
        self.solvationforces = None
        # Particle Diameters
        self.diameters = diameters

        # Combination of all diameters
        self._combinations = None
        # Pair Potential
        self._pairs = None
        # Delta Height, how mich is the particle size increased for steric repulsion (in nm)
        # The hertzian spring stiffness k_n is set individually to match the equilibrium distance
        self.delta_h = delta_h
        # Humidity
        self.humidity = humidity
        # Value of equilibrium distance between core particles in nm
        self.equi_dist = equi_dist
        # Value of normal damping divided by particle mass in kg/m^2*s
        self._gamma_n = None
        self._gamma_n0 = gamma_n
        # Value of tangential damping divided by particle mass in kg/m^2*s
        self._gamma_t = None
        self._gamma_t0 = gamma_t
        # Shear strength for sliding friction in N/m^2 as used in equation F_t = tau_s * A
        self._tau_s = None
        self._tau_s0 = tau_s
        # Shear strength for rolling friction in N/m^2 as used in equation F_r = tau_r * A
        self._tau_r0 = tau_r
        self._tau_r = None
        # Array with normal elasic coefficients
        self._kn = None
        self._etar = None
        # Waterlayer around particles
        self._waterlayer = None
        # Parameters for solvation forces
        self.sigma_solv = sigma_solv
        self.F_solv_b = F_solv_b
        self.F_solv_m = F_solv_m
        self.iswall = False
        print(self.diameters)

    @property
    def combinations(self):
        """
            Calculate the combinations of all diameter and write them into a list
            :return combinations:
        """
        if self._combinations is None:
            self._combinations = np.array([item for item in (combinations_with_replacement(self.diameters, 2))])
        return self._combinations

    @property
    def pairs(self):
        """
            Generate a list for the pair coefficients
        """
        if self._pairs is None:
            self._pairs = np.arange(1, len(self.diameters) + 1, 1)
            self._pairs = [item for item in (combinations_with_replacement(self._pairs, 2))]
        return self._pairs

    @property
    def tau_s(self):
        """
        The sliding friction. If there is a stl wall involved the value will be substituted against a
        large value

        """
        if self._tau_s is None:
            self._tau_s = np.zeros((len(self.diameters), len(self.diameters)))
            for d1 in range(len(self.diameters)):
                for d2 in range(len(self.diameters)):
                    if self.iswall and (d1 >= len(self.diameters) - 2 or
                                        d2 >= len(self.diameters) - 2):
                        self._tau_s[d1, d2] = self.slidingfrictionwall
                    else:
                        self._tau_s[d1, d2] = self._tau_s0
        return self._tau_s

    @property
    def coefficientRestitution(self):
        if self._coefficientrestitution is None:
            self._coefficientrestitution = np.zeros((len(self.diameters), len(self.diameters)))
            for d1 in range(len(self.diameters)):
                for d2 in range(len(self.diameters)):
                    self._coefficientrestitution[d1, d2] = self._coefficientrestitution0
        return self._coefficientrestitution

    @property
    def tau_r(self):
        """
        The rolling friction. If there is a stl wall involved the value will be substituted against a
        large value

        """
        if self._tau_r is None:
            self._tau_r = np.zeros((len(self.diameters), len(self.diameters)))
            for d1 in range(len(self.diameters)):
                for d2 in range(len(self.diameters)):
                    if self.iswall and (d1 >= len(self.diameters) - 2 or
                                        d2 >= len(self.diameters) - 2):
                        self._tau_r[d1, d2] = self.rollingfrictionwall
                    else:
                        self._tau_r[d1, d2] = self._tau_r0
        return self._tau_r

    @property
    def gamma_n(self):
        """
        The normal damping. If there is a stl wall involved the value will be substituted against a
        large value

        """
        if self._gamma_n is None:
            self._gamma_n = np.zeros((len(self.diameters), len(self.diameters)))
            for d1 in range(len(self.diameters)):
                for d2 in range(len(self.diameters)):
                    self._gamma_n[d1, d2] = self._gamma_n0
        return self._gamma_n

    @property
    def gamma_t(self):
        """
        The tangential damping. If there is a stl wall involved the value will be substituted against a
        large value

        """
        if self._gamma_t is None:
            self._gamma_t = np.zeros((len(self.diameters), len(self.diameters)))
            for d1 in range(len(self.diameters)):
                for d2 in range(len(self.diameters)):
                    self._gamma_t[d1, d2] = self._gamma_t0
        return self._gamma_t

    @property
    def kn(self):
        """
            Contact Stiffness in normal direction

        """
        if self._kn is None:
            self._kn = np.zeros((len(self.diameters), len(self.diameters)))
        return self._kn

    @kn.setter
    def kn(self, value):
        self._kn = value

    @property
    def waterlayer(self):
        """
        Calculate the correct water layer thickness size by linear interpolation
        from humidity vs coverage data
        return: thickness of the waterlayer
        """

        if self._waterlayer is None:
            # ---- DATA: humidity vs. coverage (10nm particles) ---- #
            data_humidity = np.array([0.0, 0.123667, 0.247692, 0.422529, 0.629865, 0.793000, 0.897458, 0.951875, 1.0])  # x
            data_waterlayer = np.array([0.262533, 0.262533, 0.330526, 0.409700, 0.486277, 0.574083, 0.697542, 0.809125, 0.809125])  # f(x)
            for lin, value in enumerate(data_humidity):
                if value > self.humidity:
                    # Linear interpolation
                    x0 = data_humidity[lin - 1]
                    x1 = data_humidity[lin]
                    f0 = data_waterlayer[lin - 1]
                    f1 = data_waterlayer[lin]
                    break
            self._waterlayer = f0 + (f1 - f0) / (x1 - x0) * (self.humidity - x0)
            print('Thickness of waterlayer: %.2f nm' % self._waterlayer)
        return self._waterlayer

    def _calculate_diameters(self, d1, d2):
        """
        Calculate the diameter of the particles
        :param d1: diameter first particle
        :param d2: diameter of second particle
        :return :
            mean diameter
            effective diemeters
            core diameter1
            core diameter2
            effective corediameter
            mean core diameter
        """
        dm = round((d1 + d2) / 2, 3)
        de = round((2 * d1 * d2) / (d1 + d2), 3)
        dc1 = round(d1 - self.equi_dist, 3)
        dc2 = round(d2 - self.equi_dist, 3)
        dce = round((2 * dc1 * dc2) / (dc1 + dc2), 3)
        dcm = round((dc1 + dc2) / 2, 3)
        print('\nMean diameter: %.2f nm' % dm)
        print('Core diameter 1: %.2f nm (without chemisorbed waterlayer of 0.3 nm thickness)' % dc1)
        print('Core diameter 2: %.2f nm (without chemisorbed waterlayer of 0.3 nm thickness)' % dc2)
        return dm, de, dc1, dc2, dce, dcm

    def _calc_surfacetension(self, dc1, dc2):
        """
        Calculate the surface tension. The surface tension of particles is chosen
        according to the next available particle size
        :param dc1: core diameter of the first particle
        :param dc2: core diameter of the second paticle
        :return surfacetension:
        """
        # ---- DATA: humidity vs. surface tension  ---- #
        # 4nm particles
        data_humidity4nm = np.array([0.0, 0.184389, 0.402225, 0.577136, 0.725398, 0.897598, 1.028767, 1.095049])
        data_gamma4nm = np.array([0.109052, 0.109052, 0.095479, 0.081527, 0.068580, 0.063857, 0.058307, 0.054609])
        # 6nm particles
        data_humidity6nm = np.array([0.0, 0.131611, 0.286943, 0.594320, 0.826887, 1.0])
        data_gamma6nm = np.array([0.107892, 0.107892, 0.086195, 0.065728, 0.060347, 0.0523])
        # 8nm particles
        data_humidity8nm = np.array([0.0, 0.134375, 0.261000, 0.578747, 0.796143, 1.0])
        data_gamma8nm = np.array([0.098067, 0.098067, 0.090194, 0.064652, 0.057125, 0.0523])
        # 10nm particles
        data_humidity10nm = np.array([0.0, 0.123667, 0.247692, 0.422529, 0.629865, 0.793000, 0.897458, 0.951875, 1.0])
        data_gamma10nm = np.array([0.081457, 0.081457, 0.067124, 0.062073, 0.055896, 0.053576, 0.051417, 0.0523, 0.0523])

        dsurf = min(dc1, dc2)

        if dsurf < 5:
            data_humidity = data_humidity4nm
            data_gamma = data_gamma4nm
        elif dsurf < 7:
            data_humidity = data_humidity6nm
            data_gamma = data_gamma6nm
        elif dsurf < 9:
            data_humidity = data_humidity8nm
            data_gamma = data_gamma8nm
        else:
            data_humidity = data_humidity10nm
            data_gamma = data_gamma10nm

        for lin, value, in enumerate(data_humidity):
            if value > self.humidity:
                x0 = data_humidity[lin - 1]
                x1 = data_humidity[lin]
                f0 = data_gamma[lin - 1]
                f1 = data_gamma[lin]
                break

        # Surface tension gets scaled from tip3p surfacetension to real surfacetension
        surfacetension = (f0 + (f1 - f0) / (x1 - x0) * (self.humidity - x0)) * (0.072 / 0.04729279)
        print("Surface Tension: %.4f nN/m" % surfacetension)
        return surfacetension

    def _calc_force_displacement_curve(self, dc1, dc2, dcm, surfacetension, contactangle, start, stop):
        """
            Run the capillary force calculation of all diameters
            :param dc1: core diameter of the first particle
            :param dc2: core diameter of the second particle
        """
        r1 = dc1 / 2 + self.waterlayer
        r2 = dc2 / 2 + self.waterlayer
        if start == 0:
            tablestart = 0
        else:
            tablestart = dcm + start
            if tablestart <= 0.01:
                print('\nWARNING: tablestart is too low, had to adjust from %f to 0.01' % tablestart)
                tablestart = 0.01
        tablestop = dcm + stop
        cforce = calc_cap_force(None,
                                theta_degree=contactangle,
                                r1=r1,
                                r2=r2,
                                humidity=self.humidity,
                                gamma=round(surfacetension, 3),
                                dist_start=tablestart,
                                dist_stop=tablestop,
                                )
        capforce = cforce.forces
        print("Kelvin length: ", cforce.lambda_k)
        return capforce

    def _calc_solvation_forces(self, reff, D):
        """
        Calculate the solvation force for the combined radius and the distances D
        :R: Combined particle radius
        :D: array ( list ) of particle-particle distances
        :return: Solvation force for the distances
        """
        F_solv_0 = self.F_solv_b + self.F_solv_m * reff
        F_solv = F_solv_0 * np.cos(2 * np.pi * D / self.sigma_solv) * np.exp(-D / self.sigma_solv)
        for ii in range(len(D)):
            # if the distance is too low, set the force to zero (particle-particle contact)
            if D[ii] < self.sigma_solv / 4:
                F_solv[ii] = 0
        return F_solv

    def _calc_combined_potential(self, capforces, solvforces, dcm, calctip):
        """
            Calculate a combined potential table for the diameters
        """
        # Generate a new combined potential table
        force = capforces[:, 1] + solvforces[:]
        if calctip:
            force *= self.tipforcescaling

        # Integrate for energy
        energy = np.zeros(len(force))
        y_total = 0
        for lin, value in enumerate(reversed(force)):
            # Multiply by 0.01 <= stepsize ("dx")
            y_total += value * 0.01
            energy[lin] = y_total
        energy = energy[::-1]
        data = [[capforces[:, 0][lin], energy[lin], value] for lin, value in enumerate(force)]
        data = np.asarray(data)
        return data, force

    def _adjust_value_in_potential(self, capforces, force, dm, d1, d2, kk, ll):
        """
            Extract the correct contact stiffness and shift the force-distance values according to the
            equilibrium distance
        """
        for lin, value in enumerate(capforces[:, 0]):
            if value == dm:
                displacement = self.delta_h
                newd1 = round(d1 + displacement, 2)
                newd2 = round(d2 + displacement, 2)
                newdm = round((newd1 + newd2) / 2, 2)
                newde = round((2 * newd1 * newd2) / (newd1 + newd2), 2)
                # R^* in LIGGGHTS is equal to Ri*Rj / (Ri+Rj) while it is equal to
                # 2*Pi*Pj / (Ri+Rj) in our case. Therefore m has to be calculated for LIGGGHTS:
                # R^* = 0.05*0.05 * newde
                m = -force[lin] / (np.sqrt(0.25 * newde) * pow(displacement, 3 / 2)) * 1e9
                print('Value of merged potential at %.2fnm (mean diameter): %.4fnN' % (dm, force[lin]))
                print('Value of kn to match equilibrium: %.2f' % (m))
                self.kn[kk, ll] = m
                self.kn[ll, kk] = m
                break
        return newd1, newd2, newdm, newde, m

    def _plotpotential(self, capforces, force, solvforces, newdm, newde, m, dm, potentialdata, d1, d2, showfig=False, savefig=True, savedat=False):
        """
        Plot the data and show it, save the figure or save the data to use it later
        """
        # Datapoints for Hertz curve
        x = np.arange(capforces[:, 0][0], newdm, 0.01)
        if np.around(x[-1], 2) != np.around(newdm, 2):
            x = np.append(x, newdm)
        y = m * 1e-9 * (0.25 * newde) ** 0.5 * (newdm - x)**(3 / 2)
        if y[-1] != y[-1]:
            y[-1] = 0
            print('WARNING: Generated "nan", adjusted it to 0')
        # Datapoints final potential
        finalPot = force[:len(x)] + y
        finalPot = np.concatenate((finalPot, force[len(x):]))

        plt.close()
        fig = plt.figure(figsize=(12, 6))
        ax = fig.add_subplot(111)
        ax.axhline(y=0, color='black', linestyle='--')
        ax.axvline(x=dm, color='black', linestyle='--')
        ax.plot(potentialdata[:, 0], potentialdata[:, 2], color='red', label=r'$F_{nb}$')
        ax.plot(x, y, lw=2, color='green', label=r'$F_n$')
        ax.plot(capforces[:, 0][:len(finalPot)], finalPot, color='orange', label=r'$F_n + F_{nb}$')
        ax.set_xlabel(r'$r$ / nm')
        ax.set_ylabel(r'$F$ / nN')
        ax.grid()
        ax.legend(loc=1)
        ax.set_title('Potential between particles with diameters %r and %r' % (d1, d2))
        if savefig or savedat:
            # Generate outputfolder for plots and empty it
            plotfolder = './Plots'
            if not os.path.exists(plotfolder):
                os.mkdir(plotfolder)
            if savefig:
                try:
                    plt.savefig(plotfolder + '/potential_%r_%r.pdf' % (d1, d2), transparent=True)
                except:
                    pass
            if savedat:
                np.savetxt(plotfolder + '/potential_Fn_%r_%r.dat' % (d1, d2), np.column_stack((x, y)))
                np.savetxt(plotfolder + '/potential_Fnb_%r_%r.dat' % (d1, d2), np.column_stack((potentialdata[:, 0], potentialdata[:, 2])))
                np.savetxt(plotfolder + '/potential_Fall_%r_%r.dat' % (d1, d2), np.column_stack((capforces[:, 0][:len(finalPot)], finalPot)))
                np.savetxt(plotfolder + '/potential_Fsolv_%r_%r.dat' % (d1, d2), np.column_stack((capforces[:, 0], solvforces)))
        if showfig:
            plt.show()

    def _add_data_to_output(self, d1, d2, newd1, newd2, surfacetension, potentialdata):
        """
        Append the potential data to outputfile
        """
        with open(self.outputfile, 'a') as f:
            f.write('\n#nominal diameter(s): %s %snm\n' % (str(d1), str(d2)))
            f.write('#new 1st DEM diameter: %.2fnm\n' % (newd1))
            f.write('#new 2nd DEM diameter: %.2fnm\n' % (newd2))
            f.write('#waterlayer: %.2fnm\t gamma %.4fnN/m\t \n' % (self.waterlayer, surfacetension))
            f.write('CAP_SOLV_%s_%s\n' % (str(d1), str(d2)))
            f.write('N %d\n\n' % (len(potentialdata)))
            for dat in range(len(potentialdata)):
                f.write('%i %e %e %e \n' % (dat + 1, potentialdata[:, 0][dat] * 1e-9,
                                            potentialdata[:, 1][dat] * 1e-18,
                                            potentialdata[:, 2][dat] * 1e-9))

    def _calculate_potential(self, d1, d2, kk, ll, contactangle, start, stop, showfig, savefig, savedat):
        """
        Calculate the potential between two particles
        :param d1: first particle diameter
        :param d2: second particle diameter
        """
        # Is one of the diameters the tip or wall?
        calc_tip = False
        calc_wall = False
        if d1 == 'tip':
            d1 = self.tipdiameter
            calc_tip = True
        if d2 == 'tip':
            d2 = self.tipdiameter
            calc_tip = True
        if d1 == 'wall':
            d1 = self.walldiameter
            calc_wall = True
        if d2 == 'wall':
            d2 = self.walldiameter
            calc_wall = True
        # Calculate effective and core diameters
        dm, de, dc1, dc2, dce, dcm = self._calculate_diameters(d1, d2)
        # Calculate the surface tension
        surfacetension = self._calc_surfacetension(dc1, dc2)
        # Calculate the capillary forces for the diameters
        capillaryforces = self._calc_force_displacement_curve(dc1, dc2, dcm,
                                                              surfacetension,
                                                              contactangle,
                                                              start,
                                                              stop)
        # Diameter for solvation forces
        D_solv = capillaryforces[:, 0] - dcm
        # Calculate solvation forces
        solvationforces = self._calc_solvation_forces(0.5 * dce, D_solv)
        # Calculate combined potential
        potentialdata, forces = self._calc_combined_potential(capillaryforces, solvationforces, dcm, calc_tip)
        # Adjust the data
        newd1, newd2, newdm, newde, m = self._adjust_value_in_potential(capillaryforces, forces, dm, d1, d2, kk, ll)
        self._plotpotential(capillaryforces, forces, solvationforces, newdm, newde, m, dm, potentialdata, d1, d2, showfig, savefig, savedat)
        self._add_data_to_output(d1, d2, newd1, newd2, surfacetension, potentialdata)
        return potentialdata

    def generate_potential_tables(self, contactangle=0, start=-1, stop=3.5, showfig=False, savefig=True, savedat=False):
        """
        Calculate the potentials for all diameter combinations and write
        the data into the outputfile. Optionally write the data into figures
        :params showfig: Show the figure for each potential (default: False)
        :params savefig: Save the figures in .png (default: True)
        :params savedat: Save the raw data for the figures (default: True)
        """
        # Generate outputfile
        with open(self.outputfile, 'w') as f:
            f.write("#potential table for nonbonded interactions: [#] [m] [eV] [N]\n#dh= %.2fnN/nm\n" % (self.delta_h))
        # Indices for the kn calculation
        kk = 0
        ll = 0
        for ii, item in enumerate(self.combinations):
            d1, d2 = item
            sys.stdout.write('\rCalculating CG Potential for diameters %2.1f %2.1f (%i/%i)        '
                             % (float(d1), float(d2), ii + 1, len(self.combinations)))
            potentialdata = self._calculate_potential(d1, d2, kk, ll, contactangle, start, stop, showfig, savefig, savedat)
            # Calculate indices for kn
            ll += 1
            if (ll >= len(self.diameters)):
                kk += 1
                ll = kk

    def generate_potential_parameter_file(self):
        """
        Generate a parameter file for LIGGGHTS which can be included into the inputfile.
        The parameterfile includes:
            kn, kt
            RollingFrictionCoefficient
            SlidingFrictionCoefficient
            Damping Rolling
            Damping Sliding
            initiation pair_style, pair_coeff
        """
        kn_string = "%d" % (len(self.diameters))
        kt_string = "%d" % (len(self.diameters))
        tau_r_string = "%d" % (len(self.diameters))
        tau_string = "%d" % (len(self.diameters))
        gamman_string = "%d" % (len(self.diameters))
        gammat_string = "%d" % (len(self.diameters))
        coefficientrestitution_string = "%d" % (len(self.diameters))
        for ii in range(len(self.kn[:, 0])):
            for jj in range(len(self.kn[0, :])):
                kn_string += " %s" % (self.kn[ii, jj])
                kt_string += " %s" % (self.kn[ii, jj])
                tau_string += " %s" % (self.tau_s[ii, jj])
                tau_r_string += " %s" % (self.tau_r[ii, jj])
                gamman_string += " %.1e" % (self.gamma_n[ii, jj])
                gammat_string += " %.1e" % (self.gamma_t[ii, jj])
                coefficientrestitution_string += " ${cRest}"  # " %.1f" %(self.coefficientRestitution[ii,jj])

        with open(self.parameterfile, 'w') as f:
            f.write("# parameter file for liggghts input created by using the following settings: humidity %s, equi_dist %s, delta_h %s, sigma_solv = %s, F_solv_b = %s, F_solv_m = %s\n" % (str(self.humidity), str(self.equi_dist), str(self.delta_h), str(self.sigma_solv), str(self.F_solv_b), str(self.F_solv_m)))
            f.write("""%s
%s
fix             p1 all property/global kn peratomtypepair %s
fix             p2 all property/global kt peratomtypepair %s
fix             p3 all property/global coefficientRollingFriction peratomtypepair %s
fix             p4 all property/global coefficientFriction peratomtypepair %s
fix             p5 all property/global gamman peratomtypepair %s
fix             p6 all property/global gammat peratomtypepair %s
pair_style      hybrid/overlay gran model hertz/stiffness/nano tangential nano cohesion bond rolling_friction nano per_molecule on table linear 450
pair_coeff      * * gran
""" % (self._delim_string('Pair Potential Parameters'), self._delim_string('Material properties required for gran pair style'), kn_string, kt_string, tau_r_string, tau_string, gamman_string, gammat_string))

            for lin, item in enumerate(self.combinations):
                f.write("pair_coeff        %i %i table %s CAP_SOLV_%s_%s\n" % (self.pairs[lin][0],
                                                                               self.pairs[lin][1],
                                                                               self.outputfile,
                                                                               str(item[0]),
                                                                               str(item[1])))

    def _delim_string(self, string, char='#'):
        length = 80
        len_part = int((length - len(string) - 2) / 2)
        if len(char * len_part + ' ' + string + ' ' + char * len_part) < length:
            return char * len_part + ' ' + string + ' ' + char * len_part + char
        else:
            return char * len_part + ' ' + string + ' ' + char * len_part

    def _delim_char(self, char='#'):
        length = 80
        return char * length
