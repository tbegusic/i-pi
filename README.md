# i-pi-eq\_neq\_spectra

This is a modification of the original i-pi code that can compute equilibrium-nonequilibrium response function associated with several types of two-time, two-dimensional spectroscopies.
The instructions below assume that you know how to use i-pi. It's documentation and tutorials can be found at https://github.com/i-pi/i-pi.
Credit for nearly all parts of this code go to the developers of i-pi. Main modifications are in the drivers/driver.f90 file and in two new files:
* drivers/pes/h2o\_dip\_pol.f90
* tools/py/eq\_neq\_spectra.py

## Installation

To get started, clone the repository, source the environment settings, and compile the driver code:
```
git clone git@github.com:tbegusic/i-pi.git

cd i-pi
source env.sh

cd drivers
make
cd ..
```

See i-pi documentation for running initial tests.

## Usage

A detailed example is provided in examples/tools/eq\_neq\_spectra. In short, the i-pi-eq\_neq\_spectra is used similar to i-pi, but takes in two inputs:
```
i-pi-eq_neq_spectra input.xml input_spec.xml > output &
```
where input.xml is the standard i-pi input and input\_spec.xml is the input that defines parameters for nonequilibrium trajectories. See example for details of the inputs.
