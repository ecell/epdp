#!/usr/bin/env python

#TODO
"""
See README for a description

# D_factor N_B N_X N

LOGLEVEL=ERROR PYTHONPATH=../.. python -O run.py 1 1 100 10

"""


from egfrd import *
from bd import *
import sys
import gfrdbase
import model

def run(outfilename, D_factor, N_B, N_X, N):
    """ 
    Function that loops over single runs of binding/rebinding
    simulations. (Plus some overhead code.)

    Arguments:
        - 

    """ #TODO
    global DX_factor, V, L, matrix_size   

    print "Printing output to: " + outfilename

    ################# Set some constants 
    radius = 2.5e-9
    sigma = radius * 2
    D = 1e-12
    D_tot = D * 2

    tau = sigma**2 / D_tot

    DX_factor = 1

    V = 1e-18 # m^3
    L = V ** (1.0/3.0) 

    kf = 0.092e-18

    DX = D * DX_factor

    matrix_size = min(max(3, int((9 * (N_X+N_B)) ** (1.0/3.0))), 60)
    print 'matrix_size=', matrix_size

    # This is the list with steps 
    T_list = [tau * .1, INF]

    # Clear the distance output file
    outfile_r = open(outfilename + '_r_-'+j+'.dat', 'a')
    outfile_r.close()

    # Open the time output file
    outfile_t = open(outfilename + '_t.dat', 'w')

    ################# Loop through the multiple (N) simulations runs
    for i in range(N):

        # Perform a single simulation run
        r_list, t_list = singlerun(T_list, D_factor, N_B, N_X)

        # Write away time data
        for t in t_list:
            outfile_t.write('%g\n' % t)
    
        # Write away distance data, write data of different 
        # measurement points to different files. (Ending up
        # with a list of distances corresponding to measure-
        # ments at times from T_list.)
        for j in range(len(r_list)):
            outfile_r = open(outfilename + '_r_-'+j+'.dat', 'a')
            outfile_r.write('%g\n' % r_list[j])
            outfile_r.close()

        # Print to screen.
        print i, r_list, t_list

        # Make sure script writes data NOW
        outfile_t.flush()

    # Close time output file
    outfile_t.close()

def singlerun(T_list, D_factor, N_B, N_X):
    """ 
    Function that performs one simulation run of binding/
    rebinding. 

    Arguments:
        -

    """ #TODO

    ################# Define the model:
    
    # Create model class
    m = model.ParticleModel(L)

    # Particle species classes
    A = model.Species('A', D, radius)
    m.add_species_type(A)
    B = model.Species('B', D, radius)
    m.add_species_type(B)
    C = model.Species('C', D, radius)
    m.add_species_type(C)
    X = model.Species('X', DX, radius)
    m.add_species_type(X)
    
    # Reaction rules
    r1 = model.create_binding_reaction_rule(A, B, C, kf)
    m.network_rules.add_reaction_rule(r1)

    r2 = model.create_unbinding_reaction_rule(C, A, B, 1e3)
    m.network_rules.add_reaction_rule(r2)

    ################# Set up simulator

    w = gfrdbase.create_world(m, matrix_size)
    nrw = gfrdbase.create_network_rules_wrapper(m)
    s = EGFRDSimulator(w, myrandom.rng, nrw)

    ################# Place particles 

    # Throw in some X particles and stirr
    if N_X != 0:
        gfrdbase.throw_in_particles(w, X, N_X)

        end_time = tau * 1
        while 1:
            s.step()
            next_time = s.get_next_time()
            if next_time > end_time:
                s.stop(end_time)
                break

    # Reset simulation, so stirring has no effect on output data
    s.reset()

    # Define positions for particles to place.
    A_pos = [0,0,0]
    B_pos = [(float(A['radius']) + float(B['radius']))+1e-23,0,0]

    # Clear an area at position A_pos and place particle A there
    """
    What happens here:
        - find particles at position A_pos within species A radius
        - delete those
        - add a number of X particles _randomly to the box_, the 
          number being equal to # removed particles.
        - if a particle is accidently placed within the "cleared"
          area, repeat the process (in very crowded box, this 
          leads to infinite loop)
        - 

    Perhaps it would be more easy to just place all the particles at
    the beginning, but give A and B initially a diffusion constant of
    0, then stirr, and then set the diffusion constant to their
    desired values..
    
    """
    print "Started clearing.."
    while 1:
        #TODO this code doesn't function properly, because you want to clear 
        # particles within a certain radius.. Now you can be removing particles 
        # that after bursting are not in the area anymore, but whose bursted 
        # domains were.

        # After bursting you should just check whether something is within the radius.

        #pp = s.get_particles_within_radius(A_pos, float(A['radius']))
        dd = s.clear_volume(A_pos, float(A['radius'])) # returns a list of 
                                                       # domains ID's that
                                                       # where bursted within
                                                       # volume
        if not dd:
            break
        for d in dd:
            for p in d.particle_id_pair:
                s.remove_particle(p)
                print "Removing particle.."
        s.throw_in_particles(X, len(pp), box1)

    place_particle(w, A, A_pos)

    # Idem for a particle B:
    while 1:
        dd = s.clear_volume(B_pos, float(B['radius'])) # returns a list of 
                                                       # domains ID's that
                                                       # where bursted within
                                                       # volume
        if not dd:
            break
        for d in dd:
            for p in d.particle_id_pair:
                s.remove_particle(p)
                print "Removing particle.."
        s.throw_in_particles(X, len(pp), box1)

    place_particle(w, B, B_pos) 
    print "Finished placing particles A & B.."

    # Place rest of B at random positions
    if N_B > 1:
        s.throw_in_particles(B, N_B-1, box1)

    r_list = []
    t_list = []
    t_last = 0

    s.step()

    next_stop = T_list[0] # In this case T_list[0] equals infinity

    i_T = 0

    #  ### Start simulating
    while 1:
        """
        What happens in this if-statement: 
        Create a list of the times particles were in bound state.

        If there was a reaction:
            - Binding: Indicated by the fact that there are no C-particles.
              In this case: record the current time (s.t) to t_last.
            - Unbinding: Indicated because this is the only other case when
              a reaction has taken place. 
              In this case: record the time binding lasted, i.e. dt between
              current time (s.t) and last reaction time (t_last).
        """
        if s.last_reaction: # aka if there was a reaction
            print s.last_reaction
            if len(s.world.get_particle_ids(C.id)) == 0:  #A,B
                print 'set t_last', s.t
                t_last = s.t  # set t_last
            else:    # C
                print 'reaction: ', s.t - t_last
                t_list.append(s.t - t_last)

        """ If it's time to log, then log distance between particles. """
        next_time = s.get_next_time()
        if next_time > next_stop:
            print 'stop', i_T, next_stop
            s.stop(next_stop)
            if len(s.world.get_particle_ids(C.id)) != 0:  #A,B
                r_list.append(0)
            # If particles are manifested as A and B, log distance
            # inbetween.
            else: # C = 0, i.e. there are particles A and B
                # r_list.append(s.distance_between_particles(A.id, B.id))
                # If there's only 1 B-particle
                particlesA = w.get_particle_ids(A.id)
                particlesB = w.get_particle_ids(B.id)
                positionA = particlesA[1].position()
                positionB = particlesB[1].position()

                particle = self.sim.world.get_particle(particle_id)[1]

                r_list.append(w.distance(positionA, positionB))

            i_T += 1
            next_stop = T_list[i_T]
        
        """ If there are no measuring moments on the list any more,
        stop """
        if next_stop == INF and len(t_list) != 0:
            print 'break', s.t
            break

        s.step()

    return r_list, t_list
    

if __name__ == '__main__':

    import os

    outfilename = 'data/rebind_' + '_'.join(sys.argv[1:4]) +\
        '_' #+ os.environ['SGE_TASK_ID']
    run(outfilename, float(sys.argv[1]), 
        int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]))




"""

Some additional notes:

    # 100 nM = 100e-9 * N_A * 100 / m^3 = 6.02e19
    # V = 1 / 6.02e19 = 1.66e-20 m^3
    # L = 2.55e-7 m

    # 1 uM = 6.02e20 / m^3
    # V = 1.66e-21 m^3
    # L = 1.18e-7

"""










