from lib.calcSolvationAndCapillaryForces import potential

# diameters in nm
diameters = [3, 4]
tablefile = 'include.potential_table'
parameterfile = 'include.potential_params'
# Target humidity
humidity = 0.5
# Damping coefficient normal contact
gamma_n = 5e6
# Damping coefficient tangential contact
gamma_t = 0
# friction coefficient sliding
tau_s = 0.65e9
# friction coefficient rolling
tau_r = 0.25e9

# Save plots of the interaction potentials
savefigures = False
# Save potential datafiles per size combination
savedata = False
# initial distance of the potential table (relative to particle diameters) in nm
start = -2
# maximum distance threshold for the potential table (relative to particle diameters) in nm
stop = 3.5


pot = potential(diameters,
                tablefile=tablefile,
                parameterfile=parameterfile,
                humidity=humidity,
                gamma_n=gamma_n,
                gamma_t=gamma_t,
                tau_s=tau_s,
                tau_r=tau_r,
                )
pot.generate_potential_tables(savefig=savefigures, savedat=savedata, start=start, stop=stop)
pot.generate_potential_parameter_file()
