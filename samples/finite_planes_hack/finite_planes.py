#!/usr/bin/python

# Small script to test transfers between planes
#
# Make sure that the egfrd system is added to your PYTHONPATH
# This means, in bash for example:
# $ export PYTHONPATH=$HOME/egfrd

# Modules
# ===============================
import sys
import os
import datetime

from egfrd import *
from visualization import vtklogger

import model
import gfrdbase
import _gfrd


# Unique seed
# ===============================
currenttime = (long(
datetime.datetime.now().month*30.5*24*3600*1e6+
datetime.datetime.now().day*24*3600*1e6+datetime.datetime.now().hour*3600*1e6+
datetime.datetime.now().minute*60*1e6+datetime.datetime.now().second*1e6+
datetime.datetime.now().microsecond))
myrandom.seed(currenttime)
print "(Seed " + str(currenttime) + ")"
 

# Constants
# ===============================
logging = False
sigma = 1e-9        # Diameter particle
D = 1e-13           # Diffusion constant
N = 20000           # Number of steps simulation will last
world_size = 1e-6   # Lengths of simulation box
NParticles = 1     # Number of particles in the simulation
k = 2e-14           # Reaction constant

# VTK Logger
# ===============================
if (logging == True):
    vtk_output_directory = 'VTK_out'
    if (os.path.exists(vtk_output_directory)):
        print '** WARNING: VTK output directory already exists, possible '+ \
                'present old VTK files might lead to errors.'


# Basic set up simulator
# ===============================
# Model
m = model.ParticleModel(world_size)

# Planes ("membranes")
plane1 = _gfrd.create_planar_surface('plane1', 
    [world_size/2,0,0], [0,0,1], [0,1,0], 
    world_size, world_size)
plane2 = _gfrd.create_planar_surface('plane2', 
    [world_size/2,0,0], [1,0,0], [0,1,0], 
    world_size, world_size)
m.add_structure(plane1)
m.add_structure(plane2)

# Species
A = model.Species('A', D, sigma, 'plane1')                                   
m.add_species_type(A)
B = model.Species('B', D, 2*sigma, 'plane2')                                   
m.add_species_type(B)

# Reaction rules
r1 = model.create_binding_reaction_rule(A, plane2, B, k)
m.network_rules.add_reaction_rule(r1)
r2 = model.create_unbinding_reaction_rule(B, plane1, A, k)
m.network_rules.add_reaction_rule(r2)


# World
w = gfrdbase.create_world(m, 3)
# Simulator
s = EGFRDSimulator(w, myrandom.rng)

# Set up logger
if (logging == True):
    l = vtklogger.VTKLogger(s, vtk_output_directory, extra_particle_step=True) 

# Throw in particles
throw_in_particles(w, A, Nparticles)
#throw_in_particles(w, B, Nparticles)

# Running 2
# ===============================
for t in range(N):

    # if logging, log and also use "try statement"
    if (logging == True): 
        # take a step (but ignore errors)
        try:
            s.step()
        except:
            print "Broken (run with logging=False for errormessage)!"
            break
        # log the step
        l.log()
    # otherwise just make steps.
    else:
        s.step()

s.stop(s.t)

if (logging == True):
    l.stop()

