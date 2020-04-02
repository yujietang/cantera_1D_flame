
import cantera as ct
import numpy as np
import sys
import csv
import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.rcParams['font.size']=15

# Simulation parameters for flameD case
p = ct.one_atm  # pressure [Pa]
width = 0.1   # m
loglevel = 1  #amount of diagnostic output (0 to 8)
MW_CH4 = 16   #mole weight of METHANE
MW_air = 29   #mole weight of air
# Tin = 294     #temperature of unburnt gas
Tmin = 240
Tmax = 400
Tstep= 20
#the mass flow rate per unit area [kg/m^2]
mmin = 0.04
mmax = 0.40
mstep = 0.01
# phi_min = 0.4 #min equivalence ratio
# phi_max = 2.05 #max equivalence ratio
# phi_step = 0.05 #equivalence step of the equivalence range [phi_min,phi_max]
tableNum = 0  #table number

plt.figure(figsize=(8,5))

#free adiabatic flame:
for Tin in np.arange(Tmin, Tmax, Tstep):
    
    #Phi=i/0.105 for CH4 of GRI3.0
    gas = ct.Solution( 'gri30.cti' )
    gas.TP = Tin, p
    gas.set_equivalence_ratio(0.9, 'CH4','O2:0.21, N2:0.78, AR:0.01')
    #reactants = 'CH4:{:}, O2:0.21, N2:0.78, AR:0.01'.format(i)  # premixed gas composition

    #Z=1/(1 + MW_air/(0.105 * Phi))

    # Set up flame object
    f = ct.FreeFlame(gas, width=width)
    f.inlet.T = Tin 
 
    f.set_refine_criteria(ratio=3, slope=0.06, curve=0.12)
    f.show_solution()
    
    # Solve with lewis number is unit
    f.transport_model = 'UnityLewis'
    f.solve(loglevel=loglevel, auto=True)
    #functions

    # def Yc():
    #     val1 = ((f.gas.Y[f.gas.species_index('CO2')])/44
    #             +(f.gas.Y[f.gas.species_index('H2O')])/18
    #             +(f.gas.Y[f.gas.species_index('H2')])/2
    #             -(f.gas.Y[f.gas.species_index('O2')])/32) 
    #     return val1

    def norm_Yc():
        # vals = np.empty(f.flame.n_points)
        val1 = ((f.gas.Y[f.gas.species_index('CO2')])/44
                +(f.gas.Y[f.gas.species_index('H2O')])/18
                +(f.gas.Y[f.gas.species_index('H2')])/2
                -(f.gas.Y[f.gas.species_index('O2')])/32)
        val1_u = ((f.Y[f.gas.species_index('CO2')][0])/44 
                +(f.Y[f.gas.species_index('H2O')][0])/18
                +(f.Y[f.gas.species_index('H2')][0])/2
                -(f.Y[f.gas.species_index('O2')][0])/32)
        val1_b = ((f.Y[f.gas.species_index('CO2')][-1])/44 
                +(f.Y[f.gas.species_index('H2O')][-1])/18
                +(f.Y[f.gas.species_index('H2')][-1])/2
                -(f.Y[f.gas.species_index('O2')][-1])/32)
        norm_val1 = (val1 - val1_u)/(val1_b - val1_u)
        # vals[i] = norm_val1
        return norm_val1
    def omega_Yc():
        """
        source term of progress variable equation
        [kmol/m^2/s]
        """
        val1 = (f.gas.net_production_rates[f.gas.species_index('CO2')]
               +f.gas.net_production_rates[f.gas.species_index('H2O')]
               +f.gas.net_production_rates[f.gas.species_index('H2')]
               -f.gas.net_production_rates[f.gas.species_index('O2')])
        return val1
    # Solve with the energy equation enabled
    f.save('CH4_adiabatic.xml', 'UnityLewis', 'solution with Le=1 transport')
    
    # write the velocity, temperature, density, and mole fractions to a CSV file
    # z = f.flame.grid
    T = f.T
    # u = f.u()
    # V = f.V()
    h = f.enthalpy_mass
    csv_file = './table_{:}.csv'.format(tableNum)
    with open(csv_file, 'w', newline = '') as outfile:
        writer = csv.writer(outfile)
        writer.writerow(list(gas.species_names) + ['OmegaYc','Yc', 'T (K)', 'h (J/kg)'])
        for n in range(f.flame.n_points):
            f.set_gas_state(n)
            # writeCSV(fcsv, [z[n], u[n], V[n], T[n], gas.density()]
            #         +list(gas.moleFractions()))
            writer.writerow(list(f.gas.Y) + [omega_Yc(), norm_Yc(), T[n], h[n]])

    f.show_solution()
    print('Le=1 flamespeed = {0:7f} m/s'.format(f.u[0]))

    # f.write_csv('./table_{:}.csv'.format(tableNum),species='Y',quiet=False)

    tableNum = tableNum + 1

    fs = 100*f.u[0]

    plt.scatter(Tin, h[0], c = 'r', marker = 'o')

#burner stablized flame:
for mdot in np.arange(mmin, mmax, mstep):
    tburner = 240.0
    gas = ct.Solution('gri30.cti')
    gas.TP = tburner, p
    gas.set_equivalence_ratio(0.9, 'CH4','O2:0.21, N2:0.78, AR:0.01')
    f = ct.BurnerFlame(gas, width=width)
    f.burner.mdot = mdot
    f.set_refine_criteria(ratio=3.0, slope=0.05, curve=0.1)
    f.transport_model = 'UnityLewis'
    f.solve(loglevel, auto=True)
    f.save('CH4_burner_flame.xml', 'Le=1', 'solution with unit Lewis number transport')
    h = f.enthalpy_mass[0] # specific enthalpy [J/kg]
    def norm_Yc():
        # vals = np.empty(f.flame.n_points)
        val1 = ((f.gas.Y[f.gas.species_index('CO2')])/44
                +(f.gas.Y[f.gas.species_index('H2O')])/18
                +(f.gas.Y[f.gas.species_index('H2')])/2
                -(f.gas.Y[f.gas.species_index('O2')])/32)
        val1_u = ((f.Y[f.gas.species_index('CO2')][0])/44 
                +(f.Y[f.gas.species_index('H2O')][0])/18
                +(f.Y[f.gas.species_index('H2')][0])/2
                -(f.Y[f.gas.species_index('O2')][0])/32)
        val1_b = ((f.Y[f.gas.species_index('CO2')][-1])/44 
                +(f.Y[f.gas.species_index('H2O')][-1])/18
                +(f.Y[f.gas.species_index('H2')][-1])/2
                -(f.Y[f.gas.species_index('O2')][-1])/32)
        norm_val1 = (val1 - val1_u)/(val1_b - val1_u)
        # vals[i] = norm_val1
        return norm_val1
    def omega_Yc():
        """
        source term of progress variable equation
        [kmol/m^2/s]
        """
        val1 = (f.gas.net_production_rates[f.gas.species_index('CO2')]
               +f.gas.net_production_rates[f.gas.species_index('H2O')]
               +f.gas.net_production_rates[f.gas.species_index('H2')]
               -f.gas.net_production_rates[f.gas.species_index('O2')])
        return val1
    # write the velocity, temperature, density, and mole fractions to a CSV file
    # z = f.flame.grid
    T = f.T
    # u = f.u()
    # V = f.V()
    h = f.enthalpy_mass
    csv_file = './table_{:}.csv'.format(tableNum)
    with open(csv_file, 'w', newline = '') as outfile:
        writer = csv.writer(outfile)
        writer.writerow(list(gas.species_names) + ['OmegaYc','Yc', 'T (K)', 'h (J/kg)'])
        for n in range(f.flame.n_points):
            f.set_gas_state(n)
            # writeCSV(fcsv, [z[n], u[n], V[n], T[n], gas.density()]
            #         +list(gas.moleFractions()))
            writer.writerow(list(f.gas.Y) + [omega_Yc(), norm_Yc(), T[n], h[n]])
    f.show_solution()

    tableNum = tableNum + 1
 
    plt.scatter(Tin, h[0], c = 'r', marker = 'o')


#plot flame speed of cantera results
plt.title('Flame speed for different value of phi')
plt.xlabel(r'T')
plt.ylabel(r'h [J/kg]')
# plt.xticks(np.linspace(0.5,2.0,5,endpoint=True))
# plt.yticks(np.linspace(0,40,9,endpoint=True))

#plot flame speed of csv file's results

#x = np.loadtxt('./comp.csv',delimiter = ',', usecols = (0,), dtype = float )
#y = np.loadtxt('./comp.csv',delimiter = ',', usecols = (1,), dtype = float )
#plt.plot(x,y)
plt.savefig('./h-T.png')

plt.show()
