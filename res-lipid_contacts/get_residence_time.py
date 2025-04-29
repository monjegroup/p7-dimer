#!/usr/bin/env python

import sys
import numpy as np
import pandas as pd
import statistics as stat
import math


# initialize dataframe
df = pd.DataFrame(columns = ["DOPC", "DPPE", "POPI", "CHOL", "DOPS"])
df_n = 0

for n in range(1,64):
        # load data file
        xy1 = np.loadtxt("dopc_" + str(n) + ".dat", delimiter=' ')
        xy2 = np.loadtxt("dppe_" + str(n) + ".dat", delimiter=' ')
        xy3 = np.loadtxt("popi_" + str(n) + ".dat", delimiter=' ')
        xy4 = np.loadtxt("chol_" + str(n) + ".dat", delimiter=' ')
        xy5 = np.loadtxt("dops_" + str(n) + ".dat", delimiter=' ')

        # sum up total lipid contacts
        cont1 = xy1[:,1]
        cont2 = xy2[:,1]
        cont3 = xy3[:,1]
        cont4 = xy4[:,1]
        cont5 = xy5[:,1]

        a = 0
        b = 0
        c = 0
        d = 0
        e = 0

        for frame in range(len(cont1)):
                if cont1[frame] >= 1:
                        a += 1
        for frame in range(len(cont2)):
                if cont2[frame] >= 1:
                        b += 1
        for frame in range(len(cont3)):
                if cont3[frame] >= 1:
                        c += 1
        for frame in range(len(cont4)):
                if cont4[frame] >= 1:
                        d += 1
        for frame in range(len(cont5)):
                if cont5[frame] >= 1:
                        e += 1

        # find residence time / percentage of times that protein is in 'contact' with each lipid
        res_time_a = (a / len(cont1)) * 100
        res_time_b = (b / len(cont2)) * 100
        res_time_c = (c / len(cont3)) * 100
        res_time_d = (d / len(cont4)) * 100
        res_time_e = (e / len(cont5)) * 100

        # append average contacts to dataframe
        df.loc[df_n] = [res_time_a, res_time_b, res_time_c, res_time_d, res_time_e]
        
        df_n += 1

# convert dataframe to numpy array
data = df.to_numpy()
print(data)

#save data
np.savetxt('pep-lipid_res_time.dat', data, delimiter='\t', fmt='%.4f', newline='\n')

