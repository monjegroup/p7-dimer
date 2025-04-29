# read coordinate and trajectory
mol new minimization_system_no_H2O.gro type gro
mol addfile prodrun_1us_sk10_no_H2O_centered.xtc type xtc first 1 last -1 step 1 waitfor all

# open file to write data
set file1 [open "helicity_p2.dat" w]

# define selection and number of frames
set sel [atomselect top "protein and (residue 63 to 125) and name CA"]
set nf [molinfo top get numframes]

# loop through each frame, collecting the helicity of the protein
for { set i 0 } { $i < $nf } { incr i } {
animate goto $i
display update ui
$sel frame $i
mol ssrecalc top
set secstruct [$sel get structure]
puts $file1 "$secstruct"
}

close $file1