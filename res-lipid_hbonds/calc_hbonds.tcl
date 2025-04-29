# Load structure and trajectory files
# edit file names (end with .gro and .xtc here) and types ('gro' and 'xtc' chosen here) as needed
mol new minimization_system_no_H2O.gro type gro
mol addfile prodrun_1us_sk10_no_H2O_centered.xtc type xtc first 1 last -1 step 1 waitfor all 

package require hbonds

##Options:

# -sel2 <atom selection> (default: none)
#-writefile <yes|no> (default: no)
# -upsel <yes|no> (update atom selections every frame? default: yes)
# -frames <begin:end> or <begin:step:end> or all or now (default: all)
#-dist <cutoff distance between donor and acceptor> (default: 3.0)
# -ang <angle cutoff> (default: 20)
#-plot <yes|no> (plot with MultiPlot, default: yes)
# -outdir <output directory> (default: current)
# -log <log filename> (default: none)
# -writefile <yes|no> (default: no)
# -outfile <dat filename> (default: hbonds.dat)
#-polar <yes|no> (consider only polar atoms (N, O, S, F)? default: no)
# -DA <D|A|both> (sel1 is the donor (D), acceptor (A), or donor and acceptor (both))
#Only valid when used with two selections, default: both)
#-type: (default: none)
#none--no detailed bonding information will be calculated
#all--hbonds in the same residue pair type are all counted
#pair--hbonds in the same residue pair type are counted once
#unique--hbonds are counted according to the donor-acceptor atom pair type
#-detailout <details output file> (default: stdout)


# lipid contacts calculations start here
for {set i 0} {$i < 63} {incr i} {
    set n [expr {$i + 1}]
    hbonds -sel1 [atomselect top "protein and residue $i"] -sel2 [atomselect top "resname DOPC"] -dist 3.2 -ang 30 -writefile yes -outfile dopc_$n.dat -plot no
    hbonds -sel1 [atomselect top "protein and residue $i"] -sel2 [atomselect top "resname DPPE"] -dist 3.2 -ang 30 -writefile yes -outfile dppe_$n.dat -plot no
    hbonds -sel1 [atomselect top "protein and residue $i"] -sel2 [atomselect top "resname POPI"] -dist 3.2 -ang 30 -writefile yes -outfile popi_$n.dat -plot no
    hbonds -sel1 [atomselect top "protein and residue $i"] -sel2 [atomselect top "resname DOPS"] -dist 3.2 -ang 30 -writefile yes -outfile dops_$n.dat -plot no
}
