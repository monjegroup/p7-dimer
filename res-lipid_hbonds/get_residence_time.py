#!/usr/bin/env python

import sys
import numpy as np
import pandas as pd
import statistics as stat
import math

#USAGE: python get_residence_time.py

# initialize dataframe
df = pd.DataFrame(columns = ["DOPC", "DPPE", "POPI", "DOPS"])
df_n = 0

for n in range(1,64):
        # load data file
        xy1 = np.loadtxt("dopc_" + str(n) + ".dat", delimiter=' ')
        xy2 = np.loadtxt("dppe_" + str(n) + ".dat", delimiter=' ')
        xy3 = np.loadtxt("popi_" + str(n) + ".dat", delimiter=' ')
        xy4 = np.loadtxt("dops_" + str(n) + ".dat", delimiter=' ')

        # sum up total lipid contacts (for last 500 ns; note there is 0.5ns saving frequency in data files)
        cont1 = xy1[-1000:,1]
        cont2 = xy2[-1000:,1]
        cont3 = xy3[-1000:,1]
        cont4 = xy4[-1000:,1]

        a = 0
        b = 0
        c = 0
        d = 0

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

        # find residence time / fraction of times that protein is in 'contact' with each lipid
        res_time_a = a / 1000
        res_time_b = b / 1000
        res_time_c = c / 1000
        res_time_d = d / 1000

        # append average contacts to dataframe
        df.loc[df_n] = [res_time_a, res_time_b, res_time_c, res_time_d]
        
        df_n += 1

# convert dataframe to numpy array
data = df.to_numpy()
print(data)

#save data
np.savetxt('hbond_res_time_last500ns.dat', data, delimiter='\t', fmt='%.4f', newline='\n')

