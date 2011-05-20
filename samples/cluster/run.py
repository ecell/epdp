
# Template file to copy and paste basic code from.

# Modules
# ===============================
import sys
#Also, set relative python path
sys.path.append('../../')
import os

import math
from egfrd import *
import model
import gfrdbase
import _gfrd
from visualization import vtklogger
#import logger
#import dumper
#from bd import *


# Constants
# ===============================
end_time = float(sys.argv[1])

sigma = 1e-5        # Diameter particle; more realistic would be e-9
D = 1e-12           # Diffusion constant
world_size = 1e-3   # Lengths of simulation box


# Functions
# ===============================
def cartesian(u1, u2, nx, ny):
#    global cartesian_x, cartesian_y, u1, u2, nx, ny
    cartesian_x = (nx * (u1[0]) + ny * (u2[0]))
    cartesian_y = (ny * (u2[1]) + nx * (u1[1]))
    return cartesian_x, cartesian_y

def in_circle(cartesian_x, cartesian_y, r_x, r_y, R):
    if (math.sqrt(math.pow((cartesian_x-r_x),2)+math.pow((cartesian_y-r_y),2)) < R):
        return True
    else:
        return False

def make_cluster(cluster_R):
    global w, A #, u1, u2, cartesian_x, cartesian_y, nx, ny

    # unit vectors make up grid to place particles on
    u1 = [sigma+1e-10,0]
    u2 = [math.cos((math.pi)/3.0)*(sigma+1e-10),math.sin((math.pi)/3.0)*(sigma+1e-10)] #TODO there is actually quite some space inbetween

    # counters 
    nx = 0
    ny = 0
    
    cartesian_x = 0
    cartesian_y = 0
    cartesian_z = 0

    # place particles on diagonal and right from it
    while (1):
        cartesian_x, cartesian_y = cartesian(u1, u2, nx, ny)
        print "placing on nx, ny: ", str(nx), ',',str(ny),\
            "(x=", str(cartesian_x),",y=", str(cartesian_y), ")"
        if in_circle(cartesian_x, cartesian_y, cluster_R, cluster_R, cluster_R):
            place_particle(w, A, [cartesian_x, cartesian_y, cartesian_z])
        ny = ny+1      
        if (cartesian_x > cluster_R*2 or cartesian_y > cluster_R*2):
            nx = nx+1
            ny = 0
            cartesian_x, cartesian_y = cartesian(u1, u2, nx, ny)
        if (cartesian_x > cluster_R*2 or cartesian_y > cluster_R*2):
            break
    # place particles left from diagonal
    print "other side"
    nx = -1
    ny = 2
    while (1):
        cartesian_x, cartesian_y = cartesian(u1, u2, nx, ny)
        print "placing on nx, ny: ", str(nx), ',',str(ny),\
            "(x=", str(cartesian_x),",y=", str(cartesian_y), ")"
        if in_circle(cartesian_x, cartesian_y, cluster_R, cluster_R, cluster_R):
            place_particle(w, A, [cartesian_x, cartesian_y, cartesian_z])
        ny = ny+1      
        if (cartesian_x > cluster_R*2 or cartesian_y > cluster_R*2):
            nx = nx-1
            ny = nx*-2
            cartesian_x, cartesian_y = cartesian(u1, u2, nx, ny)
        if (cartesian_x > cluster_R*2 or cartesian_y > cluster_R*2):
            break
  

# Basic set up simulator
# ===============================
# Model
m = model.ParticleModel(world_size)
# Species
A = model.Species('A', 0, sigma/2)                                   
m.add_species_type(A) 
B = model.Species('B', D, sigma/2)                                   
m.add_species_type(B) 
# Reaction rules
# r1 = model.create_binding_reaction_rule(A, B, X, k)
# m.network_rules.add_reaction_rule(r1)

# World
w = gfrdbase.create_world(m, 3)
# Simulator
s = EGFRDSimulator(w, myrandom.rng)

# Throw in particles
throw_in_particles(w, B, 2)

# Put in cluster
make_cluster(1e-4)

# Enable VTK Logger
# ===============================
vtk_output_directory = 'VTK_out'
if (os.path.exists(vtk_output_directory)):
    print '** WARNING: VTK output directory already exists, possible '+ \
            'present old VTK files might lead to errors.'
l = vtklogger.VTKLogger(s, vtk_output_directory, extra_particle_step=True) 

#TODO throw away:
n=0
# Running 
# ===============================
print "Setting up finished. Starting simtulation."
while ((s.get_next_time() < end_time) or (n < 1)):
    l.log()
    s.step()
    # Just making sure we have something to look at: TODO
    n = n+1

# TODO
if (n<1):
    s.stop(end_time)
else:
    s.stop(s.t)

l.stop()

print "Done, bye!"



























