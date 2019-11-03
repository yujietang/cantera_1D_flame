#modified from https://github.com/ZX114 by ZX
import numpy as np
import matplotlib.pyplot as plt
import matplotlib

phi_min = input("Enter the min equivalence ratio:\n")
phi_min = float(phi_min) #min equivalence ratio
phi_max = input("Enter the max equivalence ratio:\n")
phi_max = float(phi_max) #max equivalence ratio
phi_step = input("Enter tht equivalence step of the equivalence range [phi_min,phi_max]:\n")
phi_step = float(phi_step) #equivalence step of the equivalence range [phi_min,phi_max]
y = input("T(0) or omega(1)?\n")
filename = []
for n in np.arange(phi_min, phi_max, phi_step):
    filename.append('/home/tang/Documents/flameD/cantera/000/Phi{:.2f}.csv'.format(n))

fonts1 = 23
fonts2 = 25
linew = 3
figs = (13,9)

plt.figure(figsize=figs)
plt.rcParams['font.family'] = 'Serif'
# plt.rcParams['font.serif'] = 'Dejavu Serif'
from matplotlib import rc
plt.rcParams['mathtext.fontset'] = 'stix'

for ic,filename1 in enumerate(filename):
    data1 = np.loadtxt(filename1, delimiter=',', skiprows=1)
    data1 = np.transpose(data1)
    YAR = data1[3]
    YARO = YAR[0]
    List_fuel = data1[19]
    Yc = data1[59]
    Z1 = []
    for i in range(len(Yc)):
        Z1.append(List_fuel[0])
    
    if y == '0':
        v = [300,500,1000,1500,2000,2300]
        # norm = matplotlib.colors.Normalize(vmin=300, vmax=2300)
        sc = plt.scatter(Z1,Yc,c=data1[3],cmap='rainbow',s=2,vmin=300, vmax=2300)
        # sc = plt.scatter(Z1,Yc,c=data1[3],s=2,vmin=300, vmax=2300)
    elif y == '1':
        #v = [0.0,0.5e9,1e9,1.5e9,2e9,2.5e9]
        #norm = matplotlib.colors.Normalize(vmin=0, vmax=2.2e9)
        v = [0,5,10,15,20]
        # norm = matplotlib.colors.Normalize(vmin=0, vmax=350)
        sc = plt.scatter(Z1,Yc,c=data1[5],cmap='rainbow',s=2,vmin=0, vmax=20)

cbar = plt.colorbar(sc,ticks=v)
if y == '0':
    cbar.set_label(r'$T$ $\mathrm{(K)}$', fontsize=fonts2)
elif y == '1':
    cbar.set_label(r'$\dot \omega_{Y_c} \ \mathrm{(kg/m^3 s)}$', fontsize=fonts2)
else:
    print("Input Error")
for t in cbar.ax.get_yticklabels():
    t.set_fontsize(fonts2)
plt.xlim(0,0.155)
plt.ylim(0,1)
plt.tick_params(labelsize=fonts2)
plt.xlabel(r'$Z$ (-)',fontsize=fonts2)
plt.ylabel(r'$Y_c$ (-)',fontsize=fonts2)
plt.savefig('SLF_Z'+y+'.png',dpi=500,bbox_inches='tight')
plt.show()
