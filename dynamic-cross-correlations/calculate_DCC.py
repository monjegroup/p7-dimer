#!/usr/bin/env python

import sys
import math
import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis import align
from MDAnalysis.analysis.rms import rmsd

# USAGE on command line: python calculate_DCC.py

# This script determines the correlation between amino acid and contacting lipid movement with normalized dynamic cross correlations (nDCC) analysis.

# !!IMPORTANT NOTE!!: When the option to not normalize the DCC is chosen (i.e. 'normalized=False'), the values of DCC for each pair (i.e. each entry in the resulting matrix) is averaged over the trajectory and reported - this is what is reported in the p7 dimerization research paper.
#                     Before running this type of DCC, check that the N_occ variable is set equal to the number of frames in your trajectory. Future version will try to pass this directly into the function.


# The algorithm is adapted from https://github.com/tekpinar/correlationplus/blob/master/correlationplus/calculate.py. Hence if you use this script, Please cite:
#  Extracting Dynamical Correlations and Identifying Key Residues for Allosteric Communication in   
#  Proteins by correlationplus. Mustafa Tekpinar, Bertrand Neron, and Marc Delarue. Journal of Chemical 
#  Information and Modeling. https://doi.org/10.1021/acs.jcim.1c00742.


def writeSparseCorrData(out_file, cMatrix, \
                        selectedAtoms, \
                        Ctype: bool, 
                        symmetric: bool):
    """
        This function writes correlation data in sparse format.

        In a sparse matrix, only nonzero elements of the matrix are 
        given in 3 columns format:
        i j C_ij
        i and j are the indices of the matrix positions (or residue indices,
        not residue IDs given in PDB files).
        It returns nothing. 

    Parameters
    ----------
    out_file: string
        Correlation file to write.
    cMatrix: A numpy square matrix of floats
        Correlation matrix.
    selectedAtoms: prody object
        A list of -typically CA- atoms selected from the parsed PDB file.
    Ctype: boolean
        If Ctype=True, location indices i and j indices start from 0.
        Otherwise, it is assumed to be starting from 1.
    symmetric: boolean
        If you select it True, it will write the upper (or lower) triangle.

    Returns
    -------
    Nothing. 

    """
    data_file = open(out_file, 'w')
    n = selectedAtoms.numAtoms()
 
    if(symmetric == True):
        if (Ctype == True):
            for i in range(0, n):
                for j in range(1):
                    if(np.absolute(cMatrix[i][j]) >= 0.000001):
                        data_file.write("{0:d} {1:d} {2:.6f}\n".format(i, j, cMatrix[i][j]))
        else:
            for i in range(0, n):
                for j in range(1):
                    if(np.absolute(cMatrix[i][j]) >= 0.000001):
                        data_file.write("{0:d} {1:d} {2:.6f}\n".format((i+1), (j+1), cMatrix[i][j]))
    else:
        if (Ctype == True):
            for i in range(0, n):
                for j in range(1):
                    if(np.absolute(cMatrix[i][j]) >= 0.000001):
                        data_file.write("{0:d} {1:d} {2:.6f}\n".format(i, j, cMatrix[i][j]))
        else:
            for i in range(0, n):
                for j in range(1):
                    if(np.absolute(cMatrix[i][j]) >= 0.000001):
                        data_file.write("{0:d} {1:d} {2:.6f}\n".format((i+1), (j+1), cMatrix[i][j]))

    print("@> Finished writing correlation matrix in sparse format!")

    data_file.close()


def calcMDnDCC(topology, trajectory, lipid_name, cutoff, lipid_select, startingFrame=0, endingFrame=(-1),
               normalized=True, alignTrajectory=True,
               saveMatrix=True, out_file="normalized_DCC"):
    """
        Calculate normalized dynamical cross-correlations when a topology
        and a trajectory file is given. 

    Parameters
    ----------
    topology: string
        A GRO file.
    trajectory: string
        A trajectory file in dcd, xtc or trr format.
    lipid_name: string
        The name of the lipid for with which correlations with the protein is 
        to be evaluated.
    cutoff: int
        A distance around the protein to select neighboring lipid atoms. Unit
        = angstrom.
    lipid_select: string
        Selection that describes which lipids are to be used for the calculations.
    startingFrame: int
        You can specify this value if you want to exclude some initial frames
        from your cross-correlation calculations. Default value is 0.
    endingFrame: int
        You can specify this value if you want to calculate cross-correlation 
        calculations up to a certain ending frame. Default value is -1 and it 
        indicates the last frame in your trajectory.
    normalized: bool
        Default value is True and it means that the cross-correlation matrix
        will be normalized.
    alignTrajectory: bool
        Default value is True and it means that all frames in the trajectory 
        will be aligned to the initial frame.  
    saveMatrix: bool
        If True, cross-correlation matrix will be written to an output file. 
    out_file: string
        Output file name for the cross-correlation matrix. 
        Default value is normalized_DCC and the file extension is .dat. 

    Returns
    -------
    ccMatrix: A numpy square matrix of floats
        Cross-correlation matrix.

    """
    # Note the lipid for which correlation is being evaluated
    print("@> Selected lipid: " + lipid_name)
    # Create the universe (That sounds really fancy :)
    universe = mda.Universe(topology, trajectory)

    # Create an atomgroup from the protein alpha carbon 
    calphas = universe.select_atoms("protein and name CA")

    # Create an atomgroup from the contacting leaflet lipid phosphate atoms, based on the chosen cutoff
    patoms = universe.select_atoms(lipid_select + " and (cyzone " + str(cutoff) + " 30 -100 protein)", updating=True)

    N = calphas.n_atoms
    print(f"@> Parsed {N} Calpha atoms.")
    # Set your frame window for your trajectory that you want to analyze
    # startingFrame = 0
    if endingFrame == -1:
        endingFrame = universe.trajectory.n_frames
    skip = 1 

    # Perform Calpha alignment first
    if alignTrajectory:
        print("@> Aligning only Calpha atoms to the initial frame!")
        alignment = align.AlignTraj(universe, universe, select="protein and name CA", in_memory=False)
        alignment.run()

    Rvector1 = []
    Rvector2 = []

    # Iterate through the universe trajectory
    for timestep in universe.trajectory[startingFrame:endingFrame:skip]:
        Rvector1.append(calphas.positions.flatten())

        Pvec = patoms.positions
        if Pvec == []:
            Rvector2.append([np.nan, np.nan, np.nan])
        else:
            Rvector2.append(np.nanmean(Pvec, axis=0))

    N_Frames = len(Rvector1)

    R_average1 = np.mean(Rvector1, axis=0)
    R_average2 = np.nanmean(Rvector2, axis=0)

    print("@> Calculating cross-correlation matrix:")
    ccMatrix, N_occ = DCCmatrixCalculation(N, np.array(Rvector1), np.array(Rvector2), R_average1, R_average2)

    # Do the averaging
    if N_occ < N_Frames:
        print("@> WARNING: only " + str(N_occ) + " frames out of " + str(N_Frames) + " contained at least 1 " + lipid_name + " lipid in the selection. Keep in mind when intepreting DCC results.")
    ccMatrix = ccMatrix / float(N_occ)

    if normalized:
        cc_normalized = np.zeros((N, 1), np.double)
        for i in range(0, N):
            for j in range (1):
                cc_normalized[i][j] = ccMatrix[i][j] / ((ccMatrix[i][i] * ccMatrix[j][j]) ** 0.5)
        
        if saveMatrix:
            np.savetxt(out_file, cc_normalized, fmt='%.6f')

        return cc_normalized
    else:
        for i in range(0, N):
            for j in range(1):
                ccMatrix[i][j] = ccMatrix[i][j]
        if saveMatrix:
            np.savetxt(out_file, ccMatrix, fmt='%.6f')
        return ccMatrix
    

def DCCmatrixCalculation(N, Rvector1, Rvector2, R_average1, R_average2):
    """
        This function calculates upper triangle of dynamical cross-correlation
        matrix. 
    """
    ccMatrix = np.zeros((N, 1), np.double)
    N_occ = 1001 #IMPORTANT - edit N_occ to reflect total number of frames in your trajectory!
    for k in range(0, len(Rvector1)):
        if k % 100 == 0:
            print("@> Frame: " + str(k))
        deltaR1 = np.subtract(Rvector1[k], R_average1)
        deltaR2 = np.subtract(Rvector2[k], R_average2)
        if any(math.isnan(a) for a in deltaR2):
            N_occ = N_occ - 1
            continue
        else:
            for i in range(0, N):
                ind_3i = 3 * i
                for j in range (1):
                    ccMatrix[i][j] += (deltaR1[ind_3i]*deltaR2[j] +
                                    deltaR1[ind_3i + 1] * deltaR2[j + 1] +
                                    deltaR1[ind_3i + 2] * deltaR2[j + 2])
    return ccMatrix, N_occ



#=> COMMANDS TO CALCULATE DCC BEGIN HERE <=
# EDIT the following selection as needed to tell the program the residues that make up the leaflet that the protein is in contact with
contact_leaflet = "(resid 64:189 or resid 316:645 or resid 976:1041 or resid 1108:1161 or resid 1216:1239) and "
# important note: for these systems, the 2nd monomer residues also have resid 1-63, making the 1st lipid resid 64, not 127. This is minus 63 of those listed in gromacs index files.

# EDIT the following selection as needed to choose desired lipid species' P/O3 atoms for correlation calculation for current system. 
# To do this, simply replace the text after resname with each lipid in your system. In this example, the membrane contained 5 species: 'DOPC', 'DPPE', 'POPI', 'CHL1' and 'DOPS'.
pc_lipid_atoms = contact_leaflet + "resname DOPC and name P"
pe_lipid_atoms = contact_leaflet + "resname DPPE and name P"
pi_lipid_atoms = contact_leaflet + "resname POPI and name P"
chol_lipid_atoms = contact_leaflet + "resname CHL1 and name O3"
ps_lipid_atoms = contact_leaflet + "resname DOPS and name P"

# run DCC analysis with the commands below
DCC_dopc = calcMDnDCC("dcc_calculate.gro", "dcc_calculate_sk10.xtc", "DOPC", 20, pc_lipid_atoms, startingFrame=0, endingFrame=(-1), normalized=False, alignTrajectory=True, saveMatrix=True, out_file="nDCC_dopc.dat")

DCC_dppe = calcMDnDCC("dcc_calculate.gro", "dcc_calculate_sk10.xtc", "DPPE", 20, pe_lipid_atoms, startingFrame=0, endingFrame=(-1), normalized=False, alignTrajectory=True, saveMatrix=True, out_file="nDCC_dppe.dat")

DCC_popi = calcMDnDCC("dcc_calculate.gro", "dcc_calculate_sk10.xtc", "POPI", 20, pi_lipid_atoms, startingFrame=0, endingFrame=(-1), normalized=False, alignTrajectory=True, saveMatrix=True, out_file="nDCC_popi.dat")

DCC_chol = calcMDnDCC("dcc_calculate.gro", "dcc_calculate_sk10.xtc", "CHOL", 20, chol_lipid_atoms, startingFrame=0, endingFrame=(-1), normalized=False, alignTrajectory=True, saveMatrix=True, out_file="nDCC_chol.dat")

DCC_dops = calcMDnDCC("dcc_calculate.gro", "dcc_calculate_sk10.xtc", "DOPS", 20, ps_lipid_atoms, startingFrame=0, endingFrame=(-1), normalized=False, alignTrajectory=True, saveMatrix=True, out_file="nDCC_dops.dat")
