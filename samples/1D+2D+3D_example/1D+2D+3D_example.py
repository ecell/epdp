#!/usr/bin/python

#####################################################
##### An example featuring a closed box made up #####
##### from six planes ("membranes"), involving: #####
#####                                           #####
#####   - 3D diffusion inside the box           #####
#####   - 2D diffusion on the planes            #####
#####   - 1D transport on rods inside the box   #####
#####   - reactions in 3D                       #####
#####   - reactions in 2D                       #####
#####   - reactions (binding to a cap) in 1D    #####
#####   - "direct binding" between 3D and 2D    #####
#####                                           #####
##### by TR Sokolowski, 2015                    #####
#####################################################

##### Start with:
##### PYTHONPATH=[egfrd directory] LOGLEVEL=[ERROR|DEBUG] python ./1D+2D+3D_example.py [parameters]
##### where [parameters] = 
#####         [seed] (random generator seed)
#####         [surface_to_bulk_ratio] (ratio btw. no. of particles on membrane/in bulk)
#####         [cap_dwelltime] (1 / unbinding rate from the rod caps)
#####         [membrane_dwelltime] (1 / unbinding rate from membrane)
#####         [membrane_slowdown] (ratio bulk diff. const. / membrane single part. diff. const.)
#####         [complex_slowdown]  (ratio bulk diff. const. / membrane complex diff. const.)
#####
##### (in that order)


# #### Importing necessary files ####
# import egfrd libraries
from egfrd import *
from bd import *
import model
import gfrdbase
import _gfrd
import os
import math
import sys

# Set up the model
def set_constants():  

    global N_run, t_max, N_log, log_time, log_step, N_save, save_points, loadfile, ws
    global _SEED_, _VLOG_, _INFO_
    global place_bulk_particles, place_planes, place_plane_particles, place_rods, place_rod_particles, place_sinks, place_sphere
    global N_bulk_particles, N_rod_particles, N_plane_particles
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
    log_time = 0.1      # seconds between two log events of the VTK logger
    t_max    = 3600.0   # maximal bound for simulated time [sec]

    N_save = 5e4        # saving interval
    save_points = []    # define arbitrary save points (in steps) here; can be useful in debugging

    loadfile=''         # specify a file with a previous state here if desired

    if(log_step > 0):
        log_time = 0.0

    # The random seed
    _SEED_ = int(par1)

    # Some control flags
    _INFO_ = 1
    _VLOG_ = 1

    # Define flags
    place_bulk_particles	= 1
    place_planes                = 1
    place_plane_particles       = 1
    place_rods                  = 1
    place_rod_particles         = 1
    place_sinks                 = 0
    place_sphere		= 0

    # How many particles to place
    N_bulk_particles		= 40
    N_rod_particles		= 20
    N_plane_particles		= surface_to_bulk_no_ratio * N_bulk_particles


def setup_model():
  
    global m, structure_types, membrane, microtubule, cap, nr, sink
    global A, Am, At, Ac, B, Bm
    global L_cyl, R_cyl, R_nucl

    # Create an instance of the ParticleModel class
    m = model.ParticleModel(ws) 

    # Define the membranes as structure_types
    structure_types = []
    membrane = _gfrd.StructureType()
    membrane['name'] = 'membrane'
    m.add_structure_type(membrane)
    structure_types.append(membrane)

    # Define the rod/microtubule structure type
    microtubule = _gfrd.StructureType()
    microtubule['name'] = 'microtubule'
    m.add_structure_type(microtubule)
    structure_types.append(microtubule)

    # Define the cap and sink structure type
    cap = _gfrd.StructureType()
    cap['name'] = 'cap'
    m.add_structure_type(cap)
    structure_types.append(cap)

    sink = _gfrd.StructureType()
    sink['name'] = 'sink'
    m.add_structure_type(sink)
    structure_types.append(sink)

    nr = _gfrd.StructureType()
    nr['name'] = 'non-reactive'
    m.add_structure_type(nr)
    structure_types.append(nr)

    #############################
    ##### General parameters ####
    #############################
    # Define species parameters
    R_part  = 30.0e-9  # Particle radius [m]
    D_part  = 1e-12    # Particle bulk diff. const. [m^2/sec]
    v_drift = 0.5e-6   # [m/sec]
    L_cyl   = 0.45*ws
    R_cyl   = 25.0e-9  # cylinder radius, 25 nm = MT radius
    R_nucl  = 0.15*ws
    
    ####################################
    ##### Define the reaction rates ####
    ####################################
    # Binding to the membrane
    kb_to_membrane   = 1e-6
    ku_from_membrane = 0
    # Which of the bulk species should bind to the membrane
    bind_A           = 1
    bind_B           = 0

    # Binding to the rod
    kb_to_rod       = 1e-11
    ku_from_rod     = 0

    # Binding to the rod cap
    kb_to_cap        = 0
    ku_from_cap      = 0    
        # should be zero if the following are not, and vice versa
    # Binding to cap cluster
    kb_to_cc         = 1e-5
    ku_from_cc       = 1.0/cap_dwelltime

    # Binding to a sink on the rod
    kb_to_sink       = 0
    kb_from_sink     = 0 # not used so far

    # Reactions on the membrane
    kb_on_membrane = 1e-11
    ku_on_membrane = 1.0/membrane_dwelltime

    # Reactions between bulk and membrane ("direct binding")
    kb_to_membrane_particle = kb_on_membrane * 2.0*Pi*R_part
    ku_from_membrane_particle = ku_on_membrane
    
    ##############################
    ##### Create the species #####
    ##############################
    # Define chemical species
    # Syntax: Species(name, D, radius, [surface], [drift])
    A   = model.Species('A',   D_part, R_part)
    Am  = model.Species('Am',  D_part/membrane_slowdown, R_part, membrane)
    At  = model.Species('At',  D_part/10, R_part, microtubule, v_drift)
    Ac  = model.Species('Ac',  0, R_part, cap)
    B   = model.Species('B',   D_part, R_part)
    Bm  = model.Species('Bm',  D_part/membrane_slowdown, R_part, membrane)
    AB  = model.Species('AB',  D_part, 1.5*R_part)
    ABm = model.Species('ABm', D_part/complex_slowdown, 1.5*R_part, membrane)    
    As  = model.Species('E',   0, R_part, sink)

    # Add these species to the previously defined particle model
    m.add_species_type(A)
    m.add_species_type(Am)
    m.add_species_type(At)
    m.add_species_type(Ac)
    m.add_species_type(B)
    m.add_species_type(Bm)
    m.add_species_type(AB)
    m.add_species_type(ABm)
    m.add_species_type(As)
    
    #######################################
    ##### Create the reaction network #####
    #######################################
    existing_rules = []
    # Binding of A to the membrane
    # Reaction rule A + membrane <-> Am(membrane)
    if(kb_to_membrane > 0.0 and bind_A):
        rb_Am = model.create_binding_reaction_rule(A, membrane, Am, kb_to_membrane)
        existing_rules.append(rb_Am)
    if(ku_from_membrane > 0.0):
        ru_Am = model.create_unimolecular_reaction_rule(Am, A, ku_from_membrane)
        existing_rules.append(ru_Am)
        
    # Binding of B to the membrane
    # Reaction rule B + membrane <-> Bm
    if(kb_to_membrane > 0.0 and bind_B):
        rb_Bm = model.create_binding_reaction_rule(B, membrane, Bm, kb_to_membrane)
        existing_rules.append(rb_Bm)
    if(ku_from_membrane > 0.0):
        ru_Bm = model.create_unimolecular_reaction_rule(Bm, B, ku_from_membrane)
        existing_rules.append(ru_Bm)

    # Binding of A to the MTs
    # Reaction rule A + microtubule <-> At(microtubule)
    if(kb_to_rod > 0.0):
        rb_At = model.create_binding_reaction_rule(A, microtubule, At, kb_to_rod)
        existing_rules.append(rb_At)
    if(ku_from_rod > 0.0):    
        ru_At = model.create_unimolecular_reaction_rule(At, A, ku_from_rod)
        existing_rules.append(ru_At)

    # Binding to the microtubule cap
    # Reaction rule At(microtubule) + cap -> Ac(cap) -> A(bulk)
    if(kb_to_cap > 0.0):
        rb_Ac = model.create_binding_reaction_rule(At, cap, Ac, kb_to_cap)
        existing_rules.append(rb_Ac)
    if(ku_from_cap > 0.0):
        ru_Ac = model.create_unimolecular_reaction_rule(Ac, A, ku_from_cap)
        existing_rules.append(ru_Ac)
    # ATTENTION! Make sure that this does not repeat the clustering reactions below!
        
    # Binding to the sink on the microtubule
    # Reaction rule At(microtubule) + sink -> As(sink) -> At(microtubule)
    if(kb_to_sink > 0.0):
        rb_As = model.create_binding_reaction_rule(At, sink, As, kb_to_sink)
        existing_rules.append(rb_As)
    if(ku_from_sink > 0.0):
        ru_As = model.create_unimolecular_reaction_rule(As, At, ku_from_sink)
        existing_rules.append(ru_As)

    # Reaction between A and B both on the membrane
    # Reaction rule Am + Bm <-> ABm
    if(kb_on_membrane > 0.0):
        rb_ABm_mem = model.create_binding_reaction_rule(Am, Bm, ABm, kb_on_membrane)
        existing_rules.append(rb_ABm_mem)
    if(ku_on_membrane > 0.0):
        ru_ABm_mem = model.create_unbinding_reaction_rule(ABm, Am, Bm, ku_on_membrane)
        existing_rules.append(ru_ABm_mem)

    # "Direct" binding of A from bulk to membrane-bound Bm
    # Reaction rule A + Bm <-> ABm
    if(kb_to_membrane_particle > 0.0):
        rb_ABm_dir = model.create_binding_reaction_rule(A, Bm, ABm, kb_to_membrane_particle)
        existing_rules.append(rb_ABm_dir)
    if(ku_from_membrane_particle > 0.0):    
        rb_ABm_dir = model.create_unbinding_reaction_rule(ABm, A, Bm, ku_from_membrane_particle)
        existing_rules.append(ru_ABm_dir)

    # "Direct" binding of B from bulk to membrane-bound Am
    # Reaction rule Am + B <-> ABm
    if(kb_to_membrane_particle > 0.0):
        rb_ABm_dir_2 = model.create_binding_reaction_rule(Am, B, ABm, kb_to_membrane_particle)
        existing_rules.append(rb_ABm_dir_2)
    if(ku_from_membrane_particle > 0.0):    
        rb_ABm_dir_2 = model.create_unbinding_reaction_rule(ABm, Am, B, ku_from_membrane_particle)
        existing_rules.append(ru_ABm_dir_2)

    ### TODO: A + B <-> C in bulk        

    # Create a cascade of cap bound species and reactions between them
    # ATTENTION! Make sure that this does not repeat the cap-binding reaction above!
    ### TODO: polish this
    Ncluster_max = 50
    Ac = []
    for i in range(0,Ncluster_max):

        cluster_radius = R_part * (1.0+1.0*i/Ncluster_max)
        Ac.append( model.Species('Ac'+str(i+1), 0, cluster_radius, cap) )
        m.add_species_type(Ac[i])

    existing_rules.append( model.create_binding_reaction_rule(At, cap, Ac[0], kb_to_cap_cluster) )
    existing_rules.append( model.create_unimolecular_reaction_rule(Ac[0], A, ku_from_cap_cluster) )

    for i in range(1,Ncluster_max):
        existing_rules.append( model.create_binding_reaction_rule(At, Ac[i-1], Ac[i], kb_to_cap_cluster) )
        existing_rules.append( model.create_unbinding_reaction_rule(Ac[i], Ac[i-1], A, ku_from_cap_cluster) )
    
    # Finalize:
    # Add all the rules to the particle model
    for rule in existing_rules:
        m.add_reaction_rule(rule)


def setup_simulator():

    global w, s, bb_1, bb_2
    
    ###############################
    ##### Set up the simulator ####
    ###############################
    # Create a world with 3x3x3 cells for the matrixsize
    w = gfrdbase.create_world(m, 3)

    # Set default bounding box
    bb_1 = [0,0,0]
    bb_2 = [0,0,0]    

    # Get default structure id
    def_sid = w.get_def_structure_id()

    # Create membranes
    if place_planes :
        create_box(w, membrane, [0.5*ws, 0.5*ws, 0.5*ws], [0.95*ws, 0.4*ws, 0.4*ws],one_sided=True)
        # redefine the bounding box for bulk particle placement
        bb_1 = [0.10*ws, 0.33*ws, 0.33*ws]
        bb_2 = [0.90*ws, 0.66*ws, 0.66*ws]

    # Create microtubules/rods with caps
    if place_rods :

        yz_frac=0.40
        ry1=yz_frac
        ry2=1.0-yz_frac
        rz1=yz_frac
        rz2=1.0-yz_frac

        rod_1f_id = create_rod(w, microtubule, cap, 'rod_1_left',  [0.5*ws, ry1*ws, rz1*ws], R_cyl, [+1,0,0], L_cyl, back_cap_structure_type=nr)
        rod_1b_id = create_rod(w, microtubule, cap, 'rod_1_right', [0.5*ws, ry1*ws, rz1*ws], R_cyl, [-1,0,0], L_cyl, back_cap_structure_type=nr)
        rod_2f_id = create_rod(w, microtubule, cap, 'rod_2_left',  [0.5*ws, ry1*ws, rz2*ws], R_cyl, [+1,0,0], L_cyl, back_cap_structure_type=nr)
        rod_2b_id = create_rod(w, microtubule, cap, 'rod_2_right', [0.5*ws, ry1*ws, rz2*ws], R_cyl, [-1,0,0], L_cyl, back_cap_structure_type=nr)
        rod_3f_id = create_rod(w, microtubule, cap, 'rod_3_left',  [0.5*ws, ry2*ws, rz1*ws], R_cyl, [+1,0,0], L_cyl, back_cap_structure_type=nr)
        rod_3b_id = create_rod(w, microtubule, cap, 'rod_3_right', [0.5*ws, ry2*ws, rz1*ws], R_cyl, [-1,0,0], L_cyl, back_cap_structure_type=nr)
        rod_4f_id = create_rod(w, microtubule, cap, 'rod_4_left',  [0.5*ws, ry2*ws, rz2*ws], R_cyl, [+1,0,0], L_cyl, back_cap_structure_type=nr)
        rod_4b_id = create_rod(w, microtubule, cap, 'rod_4_right', [0.5*ws, ry2*ws, rz2*ws], R_cyl, [-1,0,0], L_cyl, back_cap_structure_type=nr)
        
        if place_sinks :
            sink_1a = model.create_disk_surface(sink.id, 'sink_1a', [0.55*ws, rx1*ws, ry1*ws], R_cyl, [1,0,0], rod_1f_id)
            #w.add_structure(sink_1a)
            sink_1b = model.create_disk_surface(sink.id, 'sink_1b', [0.65*ws, rx1*ws, ry1*ws], R_cyl, [1,0,0], rod_1f_id)
            #w.add_structure(sink_1b)
            sink_1c = model.create_disk_surface(sink.id, 'sink_1c', [0.75*ws, rx1*ws, ry1*ws], R_cyl, [1,0,0], rod_1f_id)
            #w.add_structure(sink_1c)

    if place_sphere :

        nucleus = model.create_spherical_surface(nr.id, 'nucleus', [0.5*ws, 0.5*ws, 0.5*ws], R_nucl, def_sid)
        w.add_structure(nucleus)

    # Create an eGFRD simulator for the world
    myrandom.seed(_SEED_)
    s = EGFRDSimulator(w, myrandom.rng)

    ##########################
    ##### Place particles ####
    ##########################
    if place_planes and place_plane_particles :
        throw_in_particles(w, Bm, N_plane_particles)
        
    if place_rods and place_rod_particles :
        throw_in_particles(w, At, N_rod_particles)
        
    if place_bulk_particles:
        # this one needs a bounding box to prevent bulk particles to be placed outside of the cell
        throw_in_particles(w, A, N_bulk_particles, bb_1, bb_2)
        #throw_in_particles(w, B, N_bulk_particles/2, bb_1, bb_2)        

    # throw_in_particles() for particles on disks (Ac, As) does not work correctly yet!
    # but e.g. place_particle(w, Ac[1], [1e-6, 3.8e-6, 3.8e-6]) does

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
        print '** WARNING: VTK output directory already exists, possible '+ \
                'present old VTK files might lead to errors.'
    
    # import lib and create logger instance
    from visualization import vtklogger
    vlogger = vtklogger.VTKLogger(s, vtk_output_directory, extra_particle_step=True) 
        # VTKLogger(self, sim, dir='vtkdata', buffer_size=None, show_shells=True, 
        # extra_particle_step=True, color_dict=None)


def main():

        global input_seed, surface_to_bulk_no_ratio, cap_dwelltime, membrane_dwelltime, membrane_slowdown, complex_slowdown

        # #### Read input parameters ####
        input_seed		 = int(sys.argv[1])
        surface_to_bulk_no_ratio = float(sys.argv[2])
        cap_dwelltime		 = float(sys.argv[3])
        membrane_dwelltime	 = float(sys.argv[4])
        membrane_slowdown	 = float(sys.argv[5])
        complex_slowdown	 = float(sys.argv[6])

        assert input_seed >= 0
        assert surface_to_bulk_no_ratio >= 0.0
        assert membrane_dwelltime > 0.0
        assert membrane_slowdown > 0.0
        assert complex_slowdown > 0.0

        
        # #### Set everything up: ####
        set_constants()
        setup_model()
        setup_simulator()
        if _VLOG_:
                setup_VTK()
        print '* Initialization complete, starting simulation.'
        print '* SINGLE_SHELL_FACTOR=%s, MULTI_SHELL_FACTOR=%s, MINIMAL_SEPARATION_FACTOR=%s, SAFETY=%s, TOLERANCE=%s' \
                        % (SINGLE_SHELL_FACTOR, MULTI_SHELL_FACTOR, MINIMAL_SEPARATION_FACTOR, SAFETY, TOLERANCE)
        print '* surface_to_bulk_no_ratio=%s, cap_dwelltime=%s, membrane_dwelltime=%s, membrane_slowdown=%s, complex_slowdown=%s, input_seed=%s' \
                        % (surface_to_bulk_no_ratio, cap_dwelltime, membrane_dwelltime, membrane_slowdown, complex_slowdown, input_seed)        
        print [structure.id for structure in w.structures]
        
        
        # #### Let's roll: ####
        step = 0
        tlast = 0.0
        
        if _VLOG_:
            # log initial configuration
            vlogger.log()
        
        if loadfile != '':
        
            s.load_state(loadfile)
        
        while step < N_run:
        
            step += 1
            s.step()
            #s.check() # can be expensive

            if step % (N_run/100) == 0:
                print "***** SIMULATION PROGRESS = "+str(int(step/N_run*100.0))+"%"

            if step % N_save == 0 or step in save_points:
                print "***** SAVING STATE at step "+str(int(step))
                s.save_state('state_'+str(step)+'.ini', reload=True)
        
            if _VLOG_ and step > (N_run-N_log) and \
              ( (log_step and step%log_step==0) or \
                (log_time and s.t-tlast>=log_time) ):
        
                print "***** LOGGER: logging time = "+str(s.t)
                tlast = s.t
                vlogger.log()    
                
        
        #### Quitting: ####
        print '* Done, quitting.'
        s.stop(s.t)
        #s.burst_all_domains()
        s.print_report()
        if _VLOG_:
                vlogger.stop()
                
        print '* Bye.'

        
if __name__ == "__main__":
        main()
