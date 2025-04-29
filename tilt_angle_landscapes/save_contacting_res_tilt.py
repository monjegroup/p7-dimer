#!/usr/bin/env python

import sys
import numpy as np
import pandas as pd
import statistics as stat
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib

#USAGE: python save_contacting_res_tilt.py

#-----------------------------------------#
# Part 1: collect contact frequency data
#-----------------------------------------#
df = pd.DataFrame(columns=range(63))

df_n = 0

for n in range(1,64):
	# load data file - this script requires the time series of protein residue-residue contacts to run. e.g., a mxn numpy array where m=traj frames and n=#of residues, 
	#				   such that it provides the contacts of a single p1 residue evaluated with each of the n residues of p2.
	#				   NOTE - see the sample p1r1_p2_contacts.dat file provided.

	xy1 = np.loadtxt('p1r'+str(n)+'_p2_contacts.dat', delimiter='\t')

	#nhalf = xy1[:100,:] # get freq w/ first 50 ns - USE FOR BOUND SYSTEMS!
	nhalf = xy1[-1000:,:] # get freq w/ last 500 ns  - USE FOR SURFACE SYSTEMS!

	res_list = []
	# calculate frequency of contact for each peptide 2 residue
	for res in range(63):
		a = 0
		for frame in range(len(nhalf)):
			if nhalf[frame,res] >= 1:
				a += 1
		freq = a / len(nhalf)
		res_list.append(freq)

	# append % contact frequency to dataframe
	df.loc[df_n] = np.array(res_list)
	df_n += 1

# convert dataframe to numpy array
freq_data = df.to_numpy()
freq_data = freq_data.astype(float)

np.savetxt('contacts_freq_data.dat', df.values, delimiter='\t', fmt='%.4f', newline='\n')


#------------------------------------------------------------------------------#
# Part 2: Record peptide1 residues with >= 0.9 frequency contact with peptide2
#------------------------------------------------------------------------------#
h1 = pd.DataFrame(index=range(12), columns=range(63)) # helix 1 (residues 5 to 16 = 12 total rows)
h2 = pd.DataFrame(index=range(22), columns=range(63)) # helix 2 (residues 20 to 41 = 22 total rows)
h3 = pd.DataFrame(index=range(11), columns=range(63)) # helix 3 (residues 48 to 58 = 11 total rows)

# p1 helix 1 residues that meet criteria
a = 0
for m in range(5,17):
	res_list = [m]
	for n in range(5,17):
		if freq_data[m,n] >= 0.9:
			res_list.append(n)
	h1.iloc[a,0:len(res_list)] = np.array(res_list)
	a += 1
np.savetxt('helix1_matches.dat', h1.values, delimiter='\t', fmt='%.0f', newline='\n')

# p1 helix 2 residues that meet criteria
a = 0
for m in range(20,42):
	res_list = [m]
	for n in range(20,42):
		if freq_data[m,n] >= 0.9:
			res_list.append(n)
	h2.iloc[a,0:len(res_list)] = np.array(res_list)
	a += 1
np.savetxt('helix2_matches.dat', h2.values, delimiter='\t', fmt='%.0f', newline='\n')

# p1 helix 3 residues that meet criteria
a = 0
for m in range(48,59):
	res_list = [m]
	for n in range(48,59):
		if freq_data[m,n] >= 0.9:
			res_list.append(n)
	h3.iloc[a,0:len(res_list)] = np.array(res_list)
	a += 1
np.savetxt('helix3_matches.dat', h3.values, delimiter='\t', fmt='%.0f', newline='\n')


#-----------------------------------------------------------------#
# Part 3: save tilt angles of each pair for plotting on a 2D space
#-----------------------------------------------------------------#
pep1 = np.loadtxt('pep1_tilt_angles.dat', delimiter='\t')
pep2 = np.loadtxt('pep2_tilt_angles.dat', delimiter='\t')

# initialize data lists
h1_angles = pd.DataFrame(columns=range(2))
h2_angles = pd.DataFrame(columns=range(2))
h3_angles = pd.DataFrame(columns=range(2))

# helix 1 tilt angles
p1_num = 5
i = 0
for index, row in h1.iterrows():
	for col, value in row.items():
		if col != h1.columns[0]:  # Skip the first column
			if np.isnan(value) != True:
				#print("part3 check")
				#print([p1_num-1, value-1])
				#print(np.mean(pep1[-1000:,int(p1_num-1)]))
				#print(np.mean(pep2[-1000:,value-1]))
				h1_angles.loc[i] = [np.mean(pep1[-1000:,int(p1_num-1)]), np.mean(pep2[-1000:,value-1])]
				i += 1
	p1_num += 1
np.savetxt('helix1_angle_data.dat', h1_angles.values, delimiter='\t', fmt='%.4f', newline='\n')

# helix 2 tilt angles
p1_num = 20
i = 0
for index, row in h2.iterrows():
	for col, value in row.items():
		if col != h2.columns[0]:  # Skip the first column
			if np.isnan(value) != True:
				h2_angles.loc[i] = [np.mean(pep1[-1000:,int(p1_num-1)]), np.mean(pep2[-1000:,value-1])]
				i += 1
	p1_num += 1
np.savetxt('helix2_angle_data.dat', h2_angles.values, delimiter='\t', fmt='%.4f', newline='\n')

# helix 3 tilt angles
p1_num = 48
i = 0
for index, row in h3.iterrows():
	for col, value in row.items():
		if col != h3.columns[0]:  # Skip the first column
			if np.isnan(value) != True:
				h3_angles.loc[i] = [np.mean(pep1[-1000:,int(p1_num-1)]), np.mean(pep2[-1000:,value-1])]
				i += 1
	p1_num += 1
np.savetxt('helix3_angle_data.dat', h3_angles.values, delimiter='\t', fmt='%.4f', newline='\n')

