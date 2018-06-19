import datetime
import numpy as np
import sys


class capforce:
    """
        Calculate capillary forces from iteration of circular approximation equations. Starts from a certain beta and loops to a minimum difference to P/P0.
        In the next step, the resulting beta_start from before is used as a starting point. Finds always next local minimum difference to P/P0.

        The potential calculation is done for a list of distances (dist_start to dist_stop, with interval interval).
        Several input parameters allow the customization
    """

    def __init__(self, outputfile='capforces.txt',
                 theta_degree=0,
                 r1=4,
                 r2=2,
                 humidity=0.5,
                 gamma=0.0523,
                 dist_start=6,
                 dist_stop=8,
                 interval=0.01):
        """
            Initialization of the object. It accepts several inputparameters to calculate the capillary forces.
            All paramters are parsed at this point.
            :outputfile: The outputfile is optional if a detailed description of the calculated potentials is required
            :theta_degree: The contact angle of the water at the particle surface. Is 0 for nanoparticles
            :r1: Radius of the first involved particle including chemisorbed waterlayer
            :r2: Radius of the second involved particle including chemisorbed waterlayer
            :humidity: Target humidity that determines the capillary forces
            :gamma: The surface tension of the water
            :dist_start: the lowest particle-particle distance for the potential
            :dist_stop: the largest particle-particle distance for the potential
            :inverval: the interval for the data table. Lower interval increases accuracy but increases calculation time
        """
        self.outputfile = outputfile
        self.gamma = gamma  # [nN/m] (Experimental=0.07199; TIP3P=0.0523; custom=0.07)
        self.theta_deg = theta_degree
        self.theta = np.radians(theta_degree)
        self.r1 = r1
        self.r2 = r2
        self.humidity = humidity
        self.interval = interval
        self.dist_start = dist_start
        self.dist_stop = dist_stop
        # Initial guess for beta
        self.beta_start = np.pi / 4
        # Gas constant R [N*m*mol^-1*K^-1]
        self.R = 8.314462175
        # Absolute temperature T [K]
        self.Temp = 300
        # Molar volume [m^3/mol]
        self.Vm = 0.000018047731
        self.lambda_k = 10**9 * self.gamma * self.Vm / (self.R * self.Temp)  # [nm] (Experimental=0.52  ; TIP3P=0.378; custom=0.5113; Lambda_k = 7.3044 * Gamma)
        # Effective radius
        self.reff = 2 * (self.r1 * self.r2) / (self.r1 + self.r2)
        # Initialize variable
        self._forces = None

    def _calculate_force(self, d, beta_start, interval=200, min_delta=0.0001, maxloops=2000):
        """
            calculate the capillary force by varying beta and comparing the calculated humidity to the target.
            The force is calculated for one specific distance
            :params d: distance of this inverval
            :beta_start: Initial guess for beta
            :interval: variation inverval for beta
            :min_delta: minimum deviation of targeted humidity and reached humidity
            :maxloops: maximum number of loops
            :return: distance of the particles, radii of the capillary, humidity, capillary force
        """
        # Effective particle distance (surface to surface)
        D = (d - (self.r1 + self.r2))
        # Initialize the comparison of target and calculated humidity
        delta_humid_end = 10
        delta_humid_old = 10
        # Initialize final values
        humid_end = 0
        F_end = 0
        beta_end = 0
        rm_end = 0
        lm_end = 0
        # Initialize initial beta
        beta = beta_start
        # Define interval for delta
        delta = np.pi / interval
        # Calculate approximation of F
        ii = 0
        """
            This while loop runs until the targeted humidity is reached
            or a defined number of maximum iterations is overcome
        """
        while(delta_humid_end > min_delta and ii < maxloops):
            ii += 1
            # Increase beta by interval
            beta = beta + delta
            # Calculation for 90 degrees impossible, therefore, it is changed and the rest of the
            # loop is skipped
            if(beta == np.pi / 2):
                print("beta = 90 degrees")
                beta += np.pi / 10
                continue
            # Calculate rm and lm (the geometrical parameters of the capillary)
            rm = (2 * self.reff * (1 - np.cos(beta)) + D) / (2.0 * np.cos(self.theta + beta))
            lm = self.reff * np.sin(beta) - rm * (1.0 - np.sin(self.theta + beta))
            # Again, geometrical constraints, that would cause errors
            if(abs(rm) <= 0.01 or lm <= 0.01):
                beta += np.pi / 4
                continue
            # Calculate the humidity from geometrical conditions
            humid_calc = np.exp(-self.lambda_k * (1.0 / rm - 1.0 / lm))
            # Calculate the deviation to target
            delta_humid = abs(humid_calc - self.humidity)
            # If deviation of humidity is lower, direction is right
            # The condition for the loop is reached and the final force and geometrica conditions are stored
            if(delta_humid < delta_humid_end):
                delta_humid_end = delta_humid
                humid_end = humid_calc
                F_end = -np.pi * self.gamma * self.reff * np.sin(beta) * (2.0 * np.sin(self.theta + beta) + self.reff * np.sin(beta) * (1.0 / rm - 1.0 / lm))
                beta_end = beta
                rm_end = rm
                lm_end = lm

            # If deviation of humidity is higher, change direction and decrease delta slightly
            if(delta_humid > delta_humid_old):
                delta *= -0.9
            delta_humid_old = delta_humid
        # If the calculation was not possible return an empty list
        if(delta_humid_end > 0.01 or beta_end < 0.0 or beta_end > np.pi or lm_end < 0.0):
            return [d, 0, 0, 0, 0, 0, 0, 0]
        else:
            return[d, F_end, beta_end, D, rm_end, lm_end, delta_humid_end, humid_end]

    @property
    def forces(self):
        """
            calculate the capillary force for the distance d in range dist_start,dist_stop
            Uses the method _calculate_force()
            :return F: list with all the capillary forces and geometric parameters
        """
        if self._forces is None:
            self._forces = []
            tmpBeta = self.beta_start
            loops = int((self.dist_stop - self.dist_start) / self.interval)
            for i in range(1, loops):
                dist = self.dist_start + i * self.interval
                tmpForce = self._calculate_force(dist, tmpBeta)
                # If a reasonable force was calculated, store the beta and use it for the next iteration
                if not tmpForce[-1] == 0:
                    tmpBeta = tmpForce[2]
                self._forces.append(tmpForce)
        return np.asarray(self._forces)

    def writefile(self):
        """
            Write a file containing the capillary forces and geometrical parameters
        """
        with open(self.outputfile, 'w') as write:
            write.write("#Surface_tension(Gamma) = %r; Particle_radius(R_eff) = %r; Contact_angle(Theta) = %r; Humidity(P/P0) = %r\n" % (self.gamma,
                                                                                                                                         self.reff,
                                                                                                                                         self.theta_deg,
                                                                                                                                         self.humidity))
            write.write("#%-13s %-14s %-10s %-8s %-8s %-8s %-14s %-10s\n" % ('Distance[nm]', 'Cap-Force[nN]', 'Beta[deg]', 'D[nm]', 'r[nm]', 'l[nm]', 'Delta_Humid[]', 'P/P0'))
            for lin in range(len(self.forces)):
                write.write("%-14r %-14r %-10r %-8r %-8r %-8r %-14r %-10r\n" % (round(float(self.forces[i][0]), 3),
                                                                                round(self.forces[i][1], 8),
                                                                                round(np.degrees(self.forces[i][2]), 6),
                                                                                round(self.forces[i][3], 3),
                                                                                round(self.forces[i][4], 4),
                                                                                round(self.forces[i][5], 4),
                                                                                round(self.forces[i][6], 6),
                                                                                round(self.forces[i][7], 6)))
