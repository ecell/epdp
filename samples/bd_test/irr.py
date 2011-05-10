#!/usr/bin/env python
"""
Brownian Dynamics test script

Does the same as the irreversible example, 
but allows for forcing the script to use
Brownian Dynamics (BD) only.

Script is programmed in such a way that 
choosing between BD and normal eGFRD (which 
includes BD) can be done by giving "BD" 
(Brownian Dynamics) or "GF" (Green's Functions) 
as first bash argument:

$ python irr.py [outputfile name] [time reaction is 
    allowed to run] [number of runs] [BD/GF]

e.g.:
$ python irr.py out.txt 1e-4 100 GF

"""

# import necessary libraries

import sys
# Tell Python where the main .py files are.
sys.path.append('../../')
# Alternatively one could use the bash command:
# $ export PYTHONPATH=../../

from bd import *
import gfrdbase
import model
import myrandom
import _gfrd

from egfrd import *


def run(outfilename, T, N, which_sim):
    """
    Simply runs simulation N times.

    Arguments: 
    - outfilename   : in which file to put output
    - T             : time single simulation is allowed to run
    - N             : number of times single sim is performed
    - which_sim     : whether to use eGFRD or BD simulator
    """

    outfile = open(outfilename, 'w')

    for i in range(N):
        d, t = singlerun(T, which_sim)
        outfile.write('%g\n' % d)

        print d, t
        #assert d == 0 or t == T

    outfile.close()


def singlerun(T, which_sim):
    """
    Performs a simulation that involves calculating the reaction time
    between two particles that are released next to each other.

    Input arguments:
    - T:    Time after which simulation is ended no matter what.
    - which_sim     : whether to use eGFRD or BD simulator
    """
    
    # Some constants
    sigma = 5e-9
    r0 = sigma
    D = 1e-12
    kf = 10 * sigma * D

    # ### Set up model
    m = model.ParticleModel(world_size=1e-5)

    A = model.Species('A', D/2, sigma/2)
    B = model.Species('B', D/2, sigma/2)
    C = model.Species('C', D/2, sigma/2)

    # Couple these species to the previously defined particle model
    m.add_species_type(A) # "m.add_species_type(species=A)" not
                          # allowed due boost package
    m.add_species_type(B) 
    m.add_species_type(C) 
    
    # Create reaction rules
    r1 = model.create_binding_reaction_rule(A, B, C, kf)
    m.network_rules.add_reaction_rule(r1)


    # ### Set up simulator
    w = gfrdbase.create_world(m, 3)

    # Set up the chosen simulator with an if-statement
    if (which_sim == "GF"):
        s = EGFRDSimulator(w, myrandom.rng)
                # EGFRDSimulator(self, world, 
                # rng=<_gfrd.RandomNumberGenerator 
                # object at 0xef5140>, network_rules=None)
    elif (which_sim == "BD"):
        s = BDSimulator(w, myrandom.rng,\
            _gfrd.NetworkRulesWrapper(w.model.network_rules)) 
            # BDSimulator isn't meant for user-friendly use,
            # so setting it u requires the "network_rules" 
            # (present reaction rules), which should be parsed
            # with the _gfrd.NetworkRulesWrapper(). 
    else:
        print "Incorrect argument given for simulator choice. Quitting."
        quit()
    
    # ### Place particles
    particleA = gfrdbase.place_particle(w, A, [0,0,0])
    particleB = gfrdbase.place_particle(w, B, [(float(A['radius']) + float(B['radius']))+1e-23,0,0])

    # ### Run the simulation
    while 1:
        next_time = s.t + s.dt
        if next_time > T: # stop if too much time has passed
            break
        s.step() # Perform simulation step


        if s.last_reaction: # stop if reaction has occured
            print 'reaction'
            return 0.0, s.t

    distance = s.distance_between_particles(particleA[1], particleB[1]) 
        # Give distance between particles

    return distance, s.t # return variables of interest: distance and time
        # s.t gives the current time, time that has passed in the simulator.

# go
if __name__ == '__main__':
    run(sys.argv[1], float(sys.argv[2]), int(sys.argv[3]), sys.argv[4])


# Some debugging attempts:
"""
def print_vars():
    print "############"
    for name in dir():
        myvalue = eval(name)
        print name, "is", type(name), "and is equal to ", myvalue

print_vars()
singlerun(float(sys.argv[2]), sys.argv[4])
print_vars()
singlerun(float(sys.argv[2]), sys.argv[4])
print_vars()
singlerun(float(sys.argv[2]), sys.argv[4])
"""













