# Example calculation
This quick example runs a short equilibrium trajectory and then restarts nonequilibrium trajectories from points along the equilibrium trajectory.
The nonequilibrium trajectories use the same forces but are initiated with a modified momentum, experiencing a kick at time zero in the direction of the derivative of the dipole moment.

The example uses lammps for the water force field, but any water force field can be used instead. The dipole moments, polarizabilities, and dipole moment derivatives are provided by
the added i-pi-driver option 'water\_dip\_pol'.
To run the example, just use the provided script:
```
./run.sh
```
which takes less than a minute to finish. The computed dipole moments, polarizabilities, and dipole moment derivatives are printed out (files simulation.dip\_0, simulation.pol\_0, simulation.dip\_der\_0) and can be compared to the provided bck\_\* files.

In this example, equilibrium dynamics runs for 100 steps. Nonequilibrium dynamics is initiated at every 20 steps of the equilibrium dynamics and each trajectory runs for 20 steps.
Dipoles and polarizabilities are printed every four steps. First $100/4 + 1 = 26$ entries in simulation.dip\_0 or simulation.pol\_0 are from equilibrium dynamics.
Then, $2 \times (100/20) = 10$ nonequilibrium trajectories are propagated. Each trajectory will print out $20 / 4 + 1 = 6$ dipole moments and polarizabilities into corresponding files.

You can use
```
./clean.sh
```
to delete outputs generated in the run.

Additional details about entries in input.xml and input\_spec.xml can be found in those input files. The run.sh script is also annotated.
