#!/usr/bin/python

#####################################################
##### A simple A+B <-> C reaction on a 2D plane #####
##### by TR Sokolowski, 2013                    #####
#####                                           #####
#####################################################

##### Start with:
##### PYTHONPATH=[egfrd directory] LOGLEVEL=[ERROR|DEBUG] python ./2D_example [seed] [no. of particles]


##### Importing necessary files #####
# import egfrd libraries
from egfrd import *
from bd import *
import model
import gfrdbase
import _gfrd
import os
import math


# Set up the model
def set_constants():

    global N_run, t_max, N_log, log_time, log_step, N_save, save_points, loadfile, ws
    global _SEED_, _VLOG_, _INFO_
    global N_plane_particles
    global Pi

    # Define some simulation parameters
    world_size = 10e-6  # side length of simulation box [m]
    ws = world_size     # only an abbreviation

    Pi = math.pi

    # Run and logging parameters
    N_run = 1e5         # number of simulation steps
    N_log = 1e9         # number of *last* steps it will log; if N_log>N_run, we will log N_run steps
    log_step = 0        # number of steps between two log outputs
			# if this is 0 the script will log in time intervals instead (defined below)
    log_time = 0.1     	# seconds between two log events of the VTK logger
    t_max    = 3600.0	# maximal bound for simulated time [sec]

    N_save = 5e4        # saving interval
    save_points = []	# define arbitrary save points (in steps) here; can be useful in debugging

    loadfile=''		# specify a file with a previous state here if desired

    if(log_step > 0):
        log_time = 0.0

    # The random seed
    _SEED_ = int(par1)

    # Some control flags
    _INFO_ = 1
    _VLOG_ = 1

    # How many particles to place
    N_plane_particles = int(par2)


def setup_model():

    global m, structure_types, membrane, A, B, C

    # Create an instance of the ParticleModel class
    m = model.ParticleModel(ws) 

    # Define the membrane structure_type
    structure_types = []
    membrane = _gfrd.StructureType()
    membrane['name'] = 'membrane'
    m.add_structure_type(membrane)
    structure_types.append(membrane)    

    #############################
    ##### General parameters ####
    #############################
    # Define species parameters
    R_part  = 10.0e-9   # Particle radius [m] 
    D_mem   = 0.5e-12   # Particle membrane diff. const. [um^2/s]
  
    ####################################
    ##### Define the reaction rates ####
    ####################################
    kb_intr 	= 1.0e-5 # Intrinsic rate of binding at contact;
			 # This is without the 2*pi*sigma prefactor (~1e-7 m), 
			 # which is added automatically when calculating the Green's function;
			 # Remember this is in [m/s]

    ku		= 1.0	# Unbinding rate = 1/complex lifetime [1/s]

    ##############################
    ##### Create the species #####
    ##############################
    # Define chemical species
    A = model.Species('A', D_mem, R_part, membrane)
    B = model.Species('B', D_mem, R_part, membrane)
    C = model.Species('C', D_mem/2.0, 2.0*R_part, membrane)

    # Add species to model
    for S in [A, B, C]:
	m.add_species_type(S)

    #######################################
    ##### Create the reaction network #####
    #######################################    

    rb = model.create_binding_reaction_rule(  A, B, C, kb_intr )
    ru = model.create_unbinding_reaction_rule(C, A, B, ku)

    for rule in [rb, ru]:
	m.network_rules.add_reaction_rule(rule)


def setup_simulator():

    global w, s, bb_1, bb_2

    # Create the world with 3x3x3 cells for the matrixsize
    w = gfrdbase.create_world(m, 3)

    # Get default structure id
    def_sid = w.get_def_structure_id()

    # Create membranes
    plane_z = 0.5*ws
    plane_corner = [-1.0*ws, -1.0*ws, plane_z]
    unit_x = [1, 0, 0]
    unit_y = [0, 1, 0]

    plane = model.create_planar_surface(membrane.id, 'plane', plane_corner, unit_x, unit_y, 3.0*ws, 3.0*ws, def_sid)
    w.add_structure(plane)		# Note: last argument of create_planar_surface() is the parent structure's ID
					# This is important for correct handling of unbindings from the substructures!

    # Define a bounding box for particle placement
    # If you want to restrict initial positions to a certain region change parameters here
    bb_1 = [0.01*ws, 0.01*ws, plane_z - 0.1*ws]
    bb_2 = [0.99*ws, 0.99*ws, plane_z + 0.1*ws]

    # Create a simulator for the world:
    myrandom.seed(_SEED_)
    s = EGFRDSimulator(w, myrandom.rng)

    # Randomly place the particles, using the bounding box from above
    throw_in_particles(w, A, N_plane_particles/2, bb_1, bb_2)
    throw_in_particles(w, B, N_plane_particles/2, bb_1, bb_2)

    if _INFO_:
	# Print info about created structures
	print "***** STRUCTURE TYPES *****"
	for st in structure_types:
		print ("%s \t %s" % (st.id, st['name']) )
	print "***** STRUCTURES *****"
	for struct in w.structures:
		print ("%s \t %s \t %s" % (struct.id, struct.sid, struct.name) )


def setup_VTK():

    global vlogger

    # Set up the visualization toolkit (VTK):

    # Set VTK output directory and check if output directory isn't 
    # already there
    vtk_output_directory = 'VTK_out'
    if (os.path.exists(vtk_output_directory)):
        print '***** WARNING: VTK output directory already exists, possible '+ \
                'present old VTK files might lead to errors.'
    
    # import lib and create logger instance
    from visualization import vtklogger
    vlogger = vtklogger.VTKLogger(s, vtk_output_directory, extra_particle_step=False) 
	      # VTKLogger(self, sim, dir='vtkdata', buffer_size=None, show_shells=True, 
	      # extra_particle_step=True, color_dict=None)


def main():

    global par1, par2

    # Parameters passed to script
    par1 = sys.argv[1] # used for seed
    par2 = sys.argv[2] # used for no. of plane particles

    # Make sure to convert parameters to right format when passing on!


    #### Set everything up: ####
    set_constants()
    setup_model()
    setup_simulator()
    if _VLOG_:
            setup_VTK()

    if _INFO_:
	print '* Initialization complete, starting simulation.'
	print '* SINGLE_SHELL_FACTOR=%s, MULTI_SHELL_FACTOR=%s, MINIMAL_SEPARATION_FACTOR=%s, SAFETY=%s, TOLERANCE=%s' %(SINGLE_SHELL_FACTOR, MULTI_SHELL_FACTOR, MINIMAL_SEPARATION_FACTOR, SAFETY, TOLERANCE)


    #### Let's roll: ####
    step = 0
    tlast = 0.0

    if _VLOG_:
	# log initial configuration
	vlogger.log()

    if loadfile != '':

	s.load_state(loadfile)

    while step < N_run and s.t <= t_max:

	step += 1
	s.step()
	#s.check() # still sporadically causes problems with development version

        if step % (N_run/100) == 0:
            print "***** SIMULATION PROGRESS = "+str(int(step/N_run*100.0))+"%"

        if step % N_save == 0 or step in save_points:
            print "***** SAVING STATE at step "+str(int(step))
            s.save_state('state_'+str(step)+'.ini', reload=True, delay=5.0)

        if _VLOG_ and step > (N_run-N_log) and \
          ( (log_step and step%log_step==0) or \
            (log_time and s.t-tlast>=log_time) ):

            print "***** LOGGER: logging time = "+str(s.t)
            tlast = s.t
            s.burst_all_domains()
            vlogger.log()    
            

    # #### Quitting: ####
    print '* Done, quitting.'

    s.stop(s.t)
    #s.burst_all_domains()
    s.print_report()
    if _VLOG_:
            vlogger.stop()

    print '* Bye.'


if __name__ == "__main__":

    main()
