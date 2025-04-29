#!/usr/bin/env python

import sys
import numpy as np
import pandas as pd

# USAGE on command line: percent_helicity.py [output identifier name - e.g. 'p1' or 'p2']

# define function for counting letters
def letterFrequency(filename, letter):
    
    num_hprot = []
    num_h1 = []
    num_h2 = []
    num_h3 = []

    # open file in read mode
    file = open(filename, 'r')

    text = file.readlines()

    for line in text:
    	num_hprot.append(line.count(letter)) # all protein residues
    	num_h1.append(line[8:31].count(letter)) # helix1 residues (5-16)
    	num_h2.append(line[38:81].count(letter)) # helix2 residues (20-41)
    	num_h3.append(line[94:115].count(letter)) # helix3 residues (48-58)

    df = pd.DataFrame(columns = ["prot", "h1", "h2", "h3"])
    df['prot'] = np.array(num_hprot)
    df['h1'] = np.array(num_h1)
    df['h2'] = np.array(num_h2)
    df['h3'] = np.array(num_h3)

    return df

# call the function and display the letter count
num_hres = letterFrequency('helicity_'+sys.argv[1]+'.dat', 'H')

print(num_hres)

percent_hres = pd.DataFrame(columns = ["prot", "h1", "h2", "h3"])
percent_hres['prot'] = (num_hres['prot'] / 63) * 100
percent_hres['h1'] = (num_hres['h1'] / 12) * 100
percent_hres['h2'] = (num_hres['h2'] / 22) * 100
percent_hres['h3'] = (num_hres['h3'] / 11) * 100

print(percent_hres) 

np.savetxt('num_helicity_'+sys.argv[1]+'.dat', num_hres, delimiter=',', fmt='%8.4f', newline='\n')
np.savetxt('percent_helicity_'+sys.argv[1]+'.dat', percent_hres, delimiter=',', fmt='%8.4f', newline='\n')

