#!/usr/bin/env python

import numpy as np
import MDAnalysis
import sys
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd

# USAGE on command line: tilt_pep_residue.py [gro structure file] [xtc/trr trajectory file]

# This script collect data from molecular dynamics trajectories to plot the tilt angle of a dimeric viral peptide amino acid side chains,
# defined by between a selected vector (head->tail) and a reference vector (head->tail)
# In this script, the selected vectors are based on peptide residue side chain terminating heavy atoms, and reference vectors on 
# stable residues that represent each of the 3 peptide helices' directions

# load structure and trajectory
u = MDAnalysis.Universe(sys.argv[1],sys.argv[2])

# total number of frames in trajectory
fr = int(u.trajectory.n_frames)
print("number of frames in trajectory: " + str(fr))

# Compute the angles for each molecule selected
# NOTE: THE DIRECTIONALITY OF THE VECTORS MATTERS (head/tail)

# define p7 residue dictionary containing the tail atom of each amino acid
am_acid = {
  "GLY": "name HA1"
  "ALA": "name CB",
  "VAL": "name CB",
  "LEU": "name CG",
  "ILE": "name CG1",
  "PHE": "name CZ",
  "PRO": "name CG",
  "SER": "name OG",
  "THR": "name CB",
  "TYR": "name OH",
  "LYS": "name NZ",
  "ARG": "name CZ",
  "HSD": "name CE1",
  "TRP": "name NE1",
  "ASN": "name CG"
}

# Define function to calculate the tilt angle in degrees
def theta(refa, refb, veca, vecb):
   ba = vecb - veca
   bc = refb - refa 
   theta = np.arccos(np.dot(ba,bc)/((np.linalg.norm(ba))*(np.linalg.norm(bc))))
   return np.rad2deg(theta)


# PART 1 - monomer 1 time series
df1 = pd.DataFrame(index=range(fr), columns=range(63))
for i in range(1,64):
   
   # selection statements for protein residue calculations
   res = u.select_atoms("bynum 1:994 and resid " + str(i), updating=True)
   h1 = res.select_atoms("name CA", updating=True)
   for key in am_acid.keys():
      if res.resnames[0] == key:
         t1 = res.select_atoms(am_acid[res.resnames[0]], updating=True)

   # calculate tilt angle of selected amino acids
   res_ang = []
   if i >= 5 and i <= 16: #numbers span resid of amino acids that belong to helix 1
      
      print("calculating for pep1 res " + str(i))

      refh = u.select_atoms("bynum 1:994 and resid 7 and name CA", updating=True) 
      reft = u.select_atoms("bynum 1:994 and resid 14 and name CA", updating=True)
      for ts in u.trajectory:
         ang = theta(refh.positions[0,:], reft.positions[0,:], h1.positions[0,:], t1.positions[0,:])
         res_ang.append(ang)
      
      df1.iloc[:,i-1] = res_ang
   
   elif i >= 20 and i <= 41: #numbers span resid of amino acids that belong to helix 2

      print("calculating for pep1 res " + str(i))

      refh = u.select_atoms("bynum 1:994 and resid 22 and name CA", updating=True) 
      reft = u.select_atoms("bynum 1:994 and resid 33 and name CA", updating=True)
      for ts in u.trajectory:
         ang = theta(refh.positions[0,:], reft.positions[0,:], h1.positions[0,:], t1.positions[0,:])
         res_ang.append(ang)

      df1.iloc[:,i-1] = res_ang
   
   elif i >= 48 and i <= 58: #numbers span resid of amino acids that belong to helix 3

      print("calculating for pep1 res " + str(i))

      refh = u.select_atoms("bynum 1:994 and resid 49 and name CA", updating=True) 
      reft = u.select_atoms("bynum 1:994 and resid 56 and name CA", updating=True)
      for ts in u.trajectory:
         ang = theta(refh.positions[0,:], reft.positions[0,:], h1.positions[0,:], t1.positions[0,:])
         res_ang.append(ang)
      
      df1.iloc[:,i-1] = res_ang

   np.savetxt('pep1_tilt_angles.dat', df1.values, delimiter='\t', fmt='%.4f', newline='\n')



# PART 2 - monomer 2 time series
df2 = pd.DataFrame(index=range(fr), columns=range(63))
for i in range(1,64):
   
   # selection statements for protein residue calculations
   res = u.select_atoms("bynum 995:1988 and resid " + str(i), updating=True)
   h1 = res.select_atoms("name CA", updating=True)
   for key in am_acid.keys():
      if res.resnames[0] == key:
         t1 = res.select_atoms(am_acid[res.resnames[0]], updating=True)

   # calculate tilt angle of selected amino acids
   res_ang = []
   if i >= 5 and i <= 16: #numbers span resid of amino acids that belong to helix 1
      
      print("calculating for pep2 res " + str(i))

      refh = u.select_atoms("bynum 995:1988 and resid 7 and name CA", updating=True) 
      reft = u.select_atoms("bynum 995:1988 and resid 14 and name CA", updating=True)
      for ts in u.trajectory: #UPDATE to choose trajectory slice as needed
         ang = theta(refh.positions[0,:], reft.positions[0,:], h1.positions[0,:], t1.positions[0,:])
         res_ang.append(ang)
      
      df2.iloc[:,i-1] = res_ang
   
   elif i >= 20 and i <= 41: #numbers span resid of amino acids that belong to helix 1

      print("calculating for pep2 res " + str(i))

      refh = u.select_atoms("bynum 995:1988 and resid 22 and name CA", updating=True) 
      reft = u.select_atoms("bynum 995:1988 and resid 33 and name CA", updating=True)
      for ts in u.trajectory: #UPDATE to choose trajectory slice as needed
         ang = theta(refh.positions[0,:], reft.positions[0,:], h1.positions[0,:], t1.positions[0,:])
         res_ang.append(ang)
      
      df2.iloc[:,i-1] = res_ang
   
   elif i >= 48 and i <= 58: #numbers span resid of amino acids that belong to helix 1

      print("calculating for pep2 res " + str(i))

      refh = u.select_atoms("bynum 995:1988 and resid 49 and name CA", updating=True) 
      reft = u.select_atoms("bynum 995:1988 and resid 56 and name CA", updating=True)
      for ts in u.trajectory: #UPDATE to choose trajectory slice as needed
         ang = theta(refh.positions[0,:], reft.positions[0,:], h1.positions[0,:], t1.positions[0,:])
         res_ang.append(ang)
      
      df2.iloc[:,i-1] = res_ang

   np.savetxt('pep2_tilt_angles.dat', df2.values, delimiter='\t', fmt='%.4f', newline='\n')
