#!/usr/bin/env python

import sys
import numpy as np
import pandas as pd
import math
import matplotlib.pyplot as plt
import seaborn as sns

#USAGE: python plot_pep-lip_correlation.py

# load data files
xy1 = np.loadtxt("nDCC_dopc.dat", delimiter=' ') # EDIT file path for these 5 lines as needed
xy2 = np.loadtxt("nDCC_dppe.dat", delimiter=' ')
xy3 = np.loadtxt("nDCC_popi.dat", delimiter=' ')
xy4 = np.loadtxt("nDCC_chol.dat", delimiter=' ')
xy5 = np.loadtxt("nDCC_dops.dat", delimiter=' ')

# create dataframes for each set of data you want to visualize. Here, the correlation was calculated using
#        all amino acids in both peptides, so the first and last 63 entries belong to p1 and p2 respectively.
pep1 = pd.DataFrame(columns = ["DOPC", "DPPE", "POPI", "CHOL", "DOPS"])
pep2 = pd.DataFrame(columns = ["DOPC", "DPPE", "POPI", "CHOL", "DOPS"])

# assign data to each dataframe
pep1["DOPC"] = xy1[:63]
pep2["DOPC"] = xy1[-63:]

pep1["DPPE"] = xy2[:63]
pep2["DPPE"] = xy2[-63:]

pep1["POPI"] = xy3[:63]
pep2["POPI"] = xy3[-63:]

pep1["CHOL"] = xy4[:63]
pep2["CHOL"] = xy4[-63:]

pep1["DOPS"] = xy5[:63]
pep2["DOPS"] = xy5[-63:]

# convert entries to a fraction based on maximum DCC[i,j] entry within p1 dataframe
# this loop below iterates to find maximum DCC[i,j] value for p1
num_new1 = 0
for i in range(pep1.shape[0]):
        num_new1 = num_new1
        for j in range(pep1.shape[1]):
                num = abs(pep1.iloc[i,j])
                if num > num_new1:
                        num_new1 = num
max_p1 = num_new1
pep1 = pep1.div(max_p1) # divide each entry for p1 by maximum DCC[i,j]

# convert entries to a fraction based on maximum DCC[i,j] entry within p2 dataframe
# this loop below iterates to find maximum DCC[i,j] value for p2
num_new2 = 0
for i in range(pep2.shape[0]):
        num_new2 = num_new2
        for j in range(pep2.shape[1]):
                num = abs(pep2.iloc[i,j])
                if num > num_new2:
                        num_new2 = num
max_p2 = num_new2
pep2 = pep2.div(max_p2) # divide each entry for p2 by maximum DCC[i,j]


# plot correlation maps
fig = plt.figure(figsize=(12,7))

# monomer 1 map
ax1 = plt.subplot(121)
ax1 = sns.heatmap(pep1, cmap='coolwarm', linewidth=0, vmin=-1, vmax=1)
ax1.set_xticklabels(ax1.get_xticklabels(), fontsize = '25', rotation=45)
ax1.set_yticks(np.arange(1, 66, 5))
ax1.set_yticklabels(['1', '6', '11', '16', '21', '26', '31', '36', '41', '46', '51', '56', '61'], fontsize = '25', rotation=0)
ax1.tick_params(length=6, width=2)
# use matplotlib.colorbar.Colorbar object
cbar = ax1.collections[0].colorbar
# here set the labelsize by 20
cbar.ax.tick_params(labelsize=25, length=6, width=2)
for _, spine in ax1.spines.items():
        spine.set_visible(True) 
        spine.set_linewidth(2)

# monomer 2 map
ax2 = plt.subplot(122)
ax2 = sns.heatmap(pep2, cmap='coolwarm', linewidth=0, vmin=-1, vmax=1)
ax2.set_xticklabels(ax1.get_xticklabels(), fontsize = '25', rotation=45)
ax2.set_yticks(np.arange(1, 66, 5))
ax2.set_yticklabels(['', '', '', '', '', '', '', '', '', '', '', '', ''], fontsize = '25', rotation=0)
ax2.tick_params(length=6, width=2)
# use matplotlib.colorbar.Colorbar object
cbar = ax2.collections[0].colorbar
# here set the labelsize by 20
cbar.ax.tick_params(labelsize=25, length=6, width=2)

for _, spine in ax2.spines.items():
        spine.set_visible(True) 
        spine.set_linewidth(2)

# save figure
fig.tight_layout()
# EDIT the file path as needed
plt.savefig("res_lipid_dcc.png",bbox_inches="tight",dpi=300)

