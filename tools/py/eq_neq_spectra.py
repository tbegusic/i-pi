#!/usr/bin/env python3

"""Runs equilibrium and nonequilibrium trajectories based on equilibrium-nonequilibrium method for nonlinear spectroscopy.
"""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2015 i-PI developers
# See the "licenses" directory for full license information.


import sys
import os
from copy import deepcopy
from math import ceil, floor
import xml.etree.ElementTree as et

# Check that we have the import path for this i-PI set and if not, add it.
dir_root = os.path.realpath(
    os.path.join(os.path.dirname(os.path.realpath(__file__)), '..'))
if not dir_root in sys.path:
    sys.path.insert(0, dir_root)

from ipi.utils.softexit import softexit
from ipi.engine.simulation import Simulation
import ipi.engine.outputs as eoutputs
from ipi.engine.initializer import init_chk

class EqNeqSpectra(object):

    """Class containing the details of equilibrium-nonequilibrium 2D spectra simulation.

    Attributes:
       epsilon: Magnitude of the external electric field.
       tsteps: Number of nonequilibrium dynamics steps.
    """

    def __init__(self, epsilon=0.1, tsteps=100):
        """Initializes EqNeqSpectra object.

        Args:
           epsilon: Magnitude of the external electric field.
           tsteps: Number of nonequilibrium dynamics steps.
        """
        self.epsilon = epsilon
        self.tsteps = tsteps


def read_spec_input(file_in):
    """Loads input file `file_in` and constructs an object of class EqNeqSpectra to store the relevant information.""" 
    tree = et.parse(file_in)
    root = tree.getroot()
    tsteps = int(root.find('./total_steps').text)
    epsilon = float(root.find('./epsilon').text)
    return EqNeqSpectra(epsilon, tsteps)

def modify_output_and_tsteps(sim, new_tsteps):
    """Modifies sim.outputs elements by deleting everything that is not dipole or polarizability. This is
       useful for nonequilibrium trajectories because we might want to print out more data from in the equilibrium
       trajectory calculation.

       Also stores useful information like stride and filename of checkpoint outputs, which we read back to
       initialize nonequilibrium trajectories, as well as the total number of steps of the equilibrium trajectory.
       The total number of steps is modified and set to new_tsteps.

       This procedure is invoked once after equilibrium trajectory is computed and before the first nonequilibrium 
       trajectory is launched.
    """ 
    #Have to loop over an auxiliary list of output elements, because sim.outputs is being modified inside the loop.
    outputs = list(sim.outputs[:])
    for o in outputs:
        if (type(o) is eoutputs.TrajectoryOutput):
            if o.what == "extras" and (o.extra_type == "polarizability" or o.extra_type == "dipole"):
                continue #Don't remove this element of output, we want to output dipoles and polarizabilities.
        #Store values that will help us loop over chk files.
        if (type(o) is eoutputs.CheckpointOutput):
            stride = o.stride
            chk_fn = o.filename
        sim.outputs.remove(o) #Remove everything that is not dipole or polarizability.
    #Modify tsteps and save old value, which will be used for looping over chk files.
    tsteps = sim.tsteps
    sim.tsteps = new_tsteps
    return chk_fn, stride, tsteps

def prepare_for_run(sim, file_in):
    """Reads initial q and p from a checkpoint file and resets step to zero. Invoked for each neq trajectory.""" 
    new_beads = init_chk(file_in)[0]
    sim.syslist[0].beads.q = new_beads.q
    sim.syslist[0].beads.p = new_beads.p
    sim.step = 0
    print(file_in)

def main(fn_input, fn_spec_input, options):
    """Main procedure:
       1) Reads fn_input and fn_spec_input.
       2) Runs equilibrium trajectory.
       3) Runs nonequilibrium trajectories.
       4) Closes the sockets and exits.
    """ 

    ######################################
    # My changes:
    eq_neq_spectra = read_spec_input(fn_spec_input)
    ######################################

    # construct simulation based on input file
    simulation = Simulation.load_from_xml(fn_input, request_banner=True, custom_verbosity=options.verbosity)

    # run the simulation
    simulation.run()

    ######################################
    # My changes:
    fn_chk_input, stride, tsteps_eq = modify_output_and_tsteps(simulation, eq_neq_spectra.tsteps)
    for step in range(ceil(eq_neq_spectra.tsteps/stride), floor(tsteps_eq/stride)):
        prepare_for_run(simulation, fn_chk_input + '_' + str(step))
        simulation.run()
    ######################################

    softexit.trigger(" @ SIMULATION: Exiting cleanly.")


if __name__ == '__main__':

    # TODO: Use argparse once we move to Python 2.7.

    from optparse import OptionParser

    parser = OptionParser(usage='%prog [options] <input file>',
                          description='The main i-PI executable used to run '
                                      'a simulation, given an XML input file.'
                          )

    parser.add_option('-V', '--verbosity', dest='verbosity', default=None,
                      choices=['quiet', 'low', 'medium', 'high', 'debug'],
                      help='Define the verbosity level.')

    options, args = parser.parse_args()

    # make sure that we have exactly one input file and it exists
    if len(args) == 0:
        parser.error('No input file name provided.')
    elif len(args) == 1:
        parser.error('Provide two input file names: one for i-pi, one for spectra.')
    elif len(args) > 2:
        parser.error('Provide only two input file names.')
    else:
        for fn_in in args:
            if not os.path.exists(fn_in):
                parser.error('Input file not found: {:s}'.format(fn_in))

    # Everything is ready. Go!
    main(args[0], args[1], options)
