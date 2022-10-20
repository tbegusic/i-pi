#!/bin/bash

#Set your path to lammps code.
MYLAMMPS="$HOME/Software/mylammps/build/lmp"

#Run i-pi_eq_neq_spectra. Similar to i-pi, but requires two inputs: i-pi input and spectra-related input.
i-pi-eq_neq_spectra input.xml input_spec.xml &> output &
sleep 15
#Launch lammps for forces.
$MYLAMMPS < in.lmp &> out_lammps &
#Launch i-pi-driver with option water_dip_pol for computing dipole moments, dipole derivatives, and polarizabilities.
#It also takes additional arguments with -o option. They are:
# 1. Number of beads.
# 2. Step for computing dipoles and polarizabilities (should match the stride in the input file).
# 3. Step for computing dipole derivatives (shoulds match the stride in the input file).
# 4. Number of steps for equilibrium dynamics (should match total_steps in input.xml).
# 5. Number of nonequilibrium dynamics steps (should match corr_steps in input_spec.xml).
#These arguments are sent directly to the driver so that it can avoid computing properties (e.g., dipoles and polarizabilities) when they are not needed (e.g., not printed out or used in any way).
#This feature works with arbitrary number of bead as long as the calculation is run in serial, i.e., not parallelized over beads.
i-pi-driver -u -h dipole -m water_dip_pol -o 1,4,20,100,20 &> out_dip_pol &
wait
