"""
Sandia flame D (Re=22400)
A freely-propagating, premixed CH4 flat flame with multicomponent
transport properties.
"""

import cantera as ct
import numpy as np
import sys
import csv
import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.rcParams['font.size']=15


# Simulation parameters for flameD case
p = 0.993*ct.one_atm  # pressure [Pa]
width = 0.1   # m
loglevel = 1  #amount of diagnostic output (0 to 8)
MW_CH4 = 16   #mole weight of METHANE
MW_air = 29   #mole weight of air
Tin = 294     #temperature of unburnt gas
phi_min = 0.4 #min equivalence ratio
phi_max = 2.05 #max equivalence ratio
phi_step = 0.05 #equivalence step of the equivalence range [phi_min,phi_max]
tableNum = 0  #table number

plt.figure(figsize=(8,5))

for Phi in np.arange(phi_min , phi_max , phi_step):
    
    #Phi=i/0.105 for CH4 of GRI3.0
    gas = ct.Solution( 'gri30.cti' )
    gas.TP = Tin, p
    gas.set_equivalence_ratio(Phi,'CH4','O2:0.21, N2:0.78, AR:0.01')
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
    
    # Solve with the energy equation enabled
    f.save('CH4_adiabatic.xml', 'UnityLewis', 'solution with Le=1 transport')
    f.show_solution()
    print('Le=1 flamespeed = {0:7f} m/s'.format(f.u[0]))

    f.write_csv('./table_{:}.csv'.format(tableNum),species='Y',quiet=False)

    tableNum = tableNum + 1

    fs = 100*f.u[0]

    plt.scatter(Phi, fs ,c = 'r',marker = 'o')

#plot flame speed of cantera results
plt.title('Flame speed for different value of phi')
plt.xlabel(r'Phi')
plt.ylabel(r'Flame speed [cm/s]')
plt.xticks(np.linspace(0.5,2.0,5,endpoint=True))
plt.yticks(np.linspace(0,40,9,endpoint=True))

#plot flame speed of csv file's results

#x = np.loadtxt('./comp.csv',delimiter = ',', usecols = (0,), dtype = float )
#y = np.loadtxt('./comp.csv',delimiter = ',', usecols = (1,), dtype = float )
#plt.plot(x,y)
plt.savefig('./flamespeed.png')

plt.show()
