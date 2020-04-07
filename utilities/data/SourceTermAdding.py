import numpy as np
import csv

index_min = input("Enter the min table index:\n")
index_min = int(index_min)
index_max = input("Enter the max table index:\n")
index_max = int(index_max)
index_step = input("Enter the table index step:\n")
index_step = int(index_step)
filename = []

for n in np.arange(index_min, index_max, 1):
    filename.append('./table_{:}.csv'.format(n))
    #dealing with data:
    data = np.loadtxt(filename[n], delimiter=",", skiprows=1) 
    data = np.transpose(data)
    print(filename[n])
    # print(data)
    #progress variable index position:
    Yc = data[-1]
    num = len(Yc)
    omegaYc = []
    for i in np.arange(0, num, 1):
        function = 25000.0*(Yc[i]*(1-Yc[i]))*(Yc[i]*(1-Yc[i]))
        omegaYc.append(function)
    print(omegaYc)
    omegaYc = np.array(omegaYc)
    data = np.append(data,[omegaYc],axis=0)
    
    data = np.transpose(data)
    np.savetxt(filename[n],data,delimiter=',', fmt='%f')
    # print(data)
