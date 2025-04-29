import MDAnalysis
from MDAnalysis.analysis import contacts
import sys
import numpy as np
import matplotlib.pyplot as plt

# USAGE: python residue_contact.py [system.gro] [system.xtc]

u = MDAnalysis.Universe(sys.argv[1],sys.argv[2])
print(len(u.trajectory))
 
sel_dopc = "(resname DOPC) and (name P)"
sel_dppe = "(resname DPPE) and (name P)"
sel_popi = "(resname POPI) and (name P)"
sel_chol = "(resname CHL1) and (name O3)"
sel_dops = "(resname DOPS) and (name P)"

dopc = u.select_atoms(sel_dopc)
dppe = u.select_atoms(sel_dppe)
popi = u.select_atoms(sel_popi)
chol = u.select_atoms(sel_chol)
dops = u.select_atoms(sel_dops)

for n in range(1,64):
    sel_heavy = "bynum 1:994 and resid " + str(n) + " and (name CA)" # selection of atoms in first p7 monomer
    print(sel_heavy)
    heavy = u.select_atoms(sel_heavy) 

    def contacts_within_cutoff(u, group_a, group_b, radius=14): #NOTE: may need to change the radius; here it is set to 14 angstrom
        timeseries = []
        for ts in u.trajectory:
            # calculate distances between group_a and group_b
            dist = contacts.distance_array(group_a.positions, group_b.positions)
            # determine which distances <= radius
            n_contacts = contacts.contact_matrix(dist, radius).sum()
            timeseries.append([ts.frame, n_contacts])
        return np.array(timeseries)

    hdopc = contacts_within_cutoff(u, dopc, heavy, radius=14)
    hdppe = contacts_within_cutoff(u, dppe, heavy, radius=14)
    hpopi = contacts_within_cutoff(u, popi, heavy, radius=14)
    hchol = contacts_within_cutoff(u, chol, heavy, radius=14)
    hdops = contacts_within_cutoff(u, dops, heavy, radius=14)

    file1 = open("dopc_" + str(n) + ".dat", "w")
    time = 0.0 
    for id in range(0, len(hdopc)):
        file1.write(str(round(time, 1)) + ' ' + str(hdopc[id][1]) + '\n')
        # My xtc trajectory file records a frame every 1 ns
        time += 1
    file1.close()

    file2 = open("dppe_" + str(n) + ".dat", "w")
    time = 0.0
    for id in range(0, len(hdppe)):
        file2.write(str(round(time, 1)) + ' ' + str(hdppe[id][1]) + '\n')
        # My xtc trajectory file records a frame every 1 ns
        time += 1
    file2.close()

    file3 = open("popi_" + str(n) + ".dat", "w")
    time = 0.0
    for id in range(0, len(hpopi)):
        file3.write(str(round(time, 1)) + ' ' + str(hpopi[id][1]) + '\n')
        # My xtc trajectory file records a frame every 1 ns
        time += 1
    file3.close()

    file4 = open("chol_" + str(n) + ".dat", "w")
    time = 0.0
    for id in range(0, len(hchol)):
        file4.write(str(round(time, 1)) + ' ' + str(hchol[id][1]) + '\n')
        # My xtc trajectory file records a frame every 1 ns
        time += 1
    file4.close()

    file5 = open("dops_" + str(n) + ".dat", "w")
    time = 0.0
    for id in range(0, len(hdops)):
        file5.write(str(round(time, 1)) + ' ' + str(hdops[id][1]) + '\n')
        # My xtc trajectory file records a frame every 1 ns
        time += 1
    file5.close()

    print("done with step " + str(n)) 
