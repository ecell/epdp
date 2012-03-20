#!/usr/bin/python

# Small script to test transfers between membranes
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
#myrandom.seed(currenttime)

myrandom.seed(523) # 323 for 3 part., 523 for one is good
print "(Seed " + str(currenttime) + ")"
 

# Constants
# ===============================
logging = True
#logging = True
sigma = 1e-8        # Diameter particle
DC = 1e-13           # Diffusion constant
N = 1000           # Number of steps simulation will last
world_size = 1e-6   # Lengths of simulation box
Nparticles = 1     # Number of particles in the simulation
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
m = model.ParticleModel(2.01*world_size)

# Planes ("membranes")
membrane1_type = _gfrd.StructureType()
membrane1_type['name'] = 'membrane1_type'
membrane2_type = _gfrd.StructureType()
membrane2_type['name'] = 'membrane2_type'
membrane3_type = _gfrd.StructureType()
membrane3_type['name'] = 'membrane3_type'
membrane4_type = _gfrd.StructureType()
membrane4_type['name'] = 'membrane4_type'
membrane5_type = _gfrd.StructureType()
membrane5_type['name'] = 'membrane5_type'
membrane6_type = _gfrd.StructureType()
membrane6_type['name'] = 'membrane6_type'

m.add_structure_type(membrane1_type)
m.add_structure_type(membrane2_type)
m.add_structure_type(membrane3_type)
m.add_structure_type(membrane4_type)
m.add_structure_type(membrane5_type)
m.add_structure_type(membrane6_type)

# Species
A = model.Species('A', DC, sigma, membrane1_type)
m.add_species_type(A)
B = model.Species('B', DC, sigma, membrane2_type)
m.add_species_type(B)
C = model.Species('C', DC, sigma, membrane3_type)
m.add_species_type(C)
D = model.Species('D', DC, sigma, membrane4_type)
m.add_species_type(D)
E = model.Species('E', DC, sigma, membrane5_type)
m.add_species_type(E)
F = model.Species('F', DC, sigma, membrane6_type)
m.add_species_type(F)

Statist = model.Species('Statist', 0.0, 10.0*sigma)
m.add_species_type(Statist)

# Reaction rules
rAB = model.create_binding_reaction_rule(A, membrane2_type, B, k)
m.network_rules.add_reaction_rule(rAB)
rBA = model.create_binding_reaction_rule(B, membrane1_type, A, k)
m.network_rules.add_reaction_rule(rBA)

rBC = model.create_binding_reaction_rule(B, membrane3_type, C, k)
m.network_rules.add_reaction_rule(rBC)
rCB = model.create_binding_reaction_rule(C, membrane2_type, B, k)
m.network_rules.add_reaction_rule(rCB)

rCD = model.create_binding_reaction_rule(C, membrane4_type, D, k)
m.network_rules.add_reaction_rule(rCD)
rDC = model.create_binding_reaction_rule(D, membrane3_type, C, k)
m.network_rules.add_reaction_rule(rDC)

rAD = model.create_binding_reaction_rule(A, membrane4_type, D, k)
m.network_rules.add_reaction_rule(rAD)
rDA = model.create_binding_reaction_rule(D, membrane1_type, A, k)
m.network_rules.add_reaction_rule(rDA)

rAE = model.create_binding_reaction_rule(A, membrane5_type, E, k)
m.network_rules.add_reaction_rule(rAE)
rEA = model.create_binding_reaction_rule(E, membrane1_type, A, k)
m.network_rules.add_reaction_rule(rEA)

rBE = model.create_binding_reaction_rule(B, membrane5_type, E, k)
m.network_rules.add_reaction_rule(rBE)
rEB = model.create_binding_reaction_rule(E, membrane2_type, B, k)
m.network_rules.add_reaction_rule(rEB)

rCE = model.create_binding_reaction_rule(C, membrane5_type, E, k)
m.network_rules.add_reaction_rule(rCE)
rEC = model.create_binding_reaction_rule(E, membrane3_type, C, k)
m.network_rules.add_reaction_rule(rEC)

rDE = model.create_binding_reaction_rule(D, membrane5_type, E, k)
m.network_rules.add_reaction_rule(rDE)
rED = model.create_binding_reaction_rule(E, membrane4_type, D, k)
m.network_rules.add_reaction_rule(rED)

rAF = model.create_binding_reaction_rule(A, membrane6_type, F, k)
m.network_rules.add_reaction_rule(rAF)
rFA = model.create_binding_reaction_rule(F, membrane1_type, A, k)
m.network_rules.add_reaction_rule(rFA)

rBF = model.create_binding_reaction_rule(B, membrane6_type, F, k)
m.network_rules.add_reaction_rule(rBF)
rFB = model.create_binding_reaction_rule(F, membrane2_type, B, k)
m.network_rules.add_reaction_rule(rFB)

rCF = model.create_binding_reaction_rule(C, membrane6_type, F, k)
m.network_rules.add_reaction_rule(rCF)
rFC = model.create_binding_reaction_rule(F, membrane3_type, C, k)
m.network_rules.add_reaction_rule(rFC)

rDF = model.create_binding_reaction_rule(D, membrane6_type, F, k)
m.network_rules.add_reaction_rule(rDF)
rFD = model.create_binding_reaction_rule(F, membrane4_type, D, k)
m.network_rules.add_reaction_rule(rFD)

rEF = model.create_binding_reaction_rule(E, membrane6_type, F, k)
m.network_rules.add_reaction_rule(rEF)
rFE = model.create_binding_reaction_rule(F, membrane5_type, E, k)
m.network_rules.add_reaction_rule(rFE)


# World
w = gfrdbase.create_world(m, 3)

# Create instances of membrane types
membraneA = model.create_planar_surface(membrane1_type.id, 'name1', [0,0,0], [1,0,0], [0,1,0], world_size, world_size)
w.add_structure(membraneA)
membraneB = model.create_planar_surface(membrane2_type.id, 'name2', [0,0,0], [0,0,1], [0,1,0], world_size, world_size)
w.add_structure(membraneB)
membraneC = model.create_planar_surface(membrane3_type.id, 'name3', [0,0,world_size], [1,0,0], [0,1,0], world_size, world_size)
w.add_structure(membraneC)
membraneD = model.create_planar_surface(membrane4_type.id, 'name4', [world_size,0,0], [0,0,1], [0,1,0], world_size, world_size)
w.add_structure(membraneD)
membraneE = model.create_planar_surface(membrane5_type.id, 'name5', [0,0,0], [1,0,0], [0,0,1], world_size, world_size)
w.add_structure(membraneE)
membraneF = model.create_planar_surface(membrane6_type.id, 'name6', [0,world_size,0], [1,0,0], [0,0,1], world_size, world_size)
w.add_structure(membraneF)

print("Membrane A : unit_z = " + str(membraneA.shape.unit_z))
print("Membrane B : unit_z = " + str(membraneB.shape.unit_z))
print("Membrane C : unit_z = " + str(membraneC.shape.unit_z))
print("Membrane D : unit_z = " + str(membraneD.shape.unit_z))
print("Membrane E : unit_z = " + str(membraneE.shape.unit_z))
print("Membrane F : unit_z = " + str(membraneF.shape.unit_z))


# Simulator
s = EGFRDSimulator(w, myrandom.rng)

# Set up logger
if (logging == True):
    l = vtklogger.VTKLogger(s, vtk_output_directory, extra_particle_step=True) 

# Throw in particles
place_particle(w, A, [world_size/2, world_size/2, 0])  # place one particle on membraneA
#place_particle(w, A, [world_size/4, world_size/4, 0])  # place one particle on membraneA
#place_particle(w, A, [3*world_size/4, 3*world_size/4, 0])  # place one particle on membraneA

# Statist particles in the corners of the system
# does not really work because they intersect with membranes...
#place_particle(w, Statist, [0, 0, 0])
#place_particle(w, Statist, [world_size, 0, 0])
#place_particle(w, Statist, [0, world_size, 0])
#place_particle(w, Statist, [0, 0, world_size])
#place_particle(w, Statist, [0, world_size, world_size])
#place_particle(w, Statist, [world_size, 0, world_size])
#place_particle(w, Statist, [world_size, world_size, 0])
#place_particle(w, Statist, [world_size, world_size, world_size])


#throw_in_particles(w, A, Nparticles) # THIS DOES NOT WORK BECAUSE PARTICLES INTERSECT WITH PLANES YET!
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

