#!/usr/bin/env python

'''
# D_factor N_B N_X N

LOGLEVEL=ERROR PYTHONPATH=../.. python -O run.py 1 1 100 10

'''


from egfrd import *
from bd import *
import sys
import gfrdbase
import model

def run(outfilename, D_factor, N_B, N_X, N):
    print outfilename

    radius = 2.5e-9
    sigma = radius * 2
    D = 1e-12
    D_tot = D * 2

    tau = sigma**2 / D_tot
    print 'tau=', tau

    #T_list = [tau * .1, INF]
    T_list = [INF]

    outfile_t = open(outfilename + '_t.dat', 'w')
    #outfile_r_list = [open(outfilename + '_r_-1.dat', 'w')] 

    for i in range(N):
        r_list, t_list = singlerun(T_list, D_factor, N_B, N_X)

        for t in t_list:
            outfile_t.write('%g\n' % t)

        #for j in range(len(r_list)):
        #    outfile_r_list[j].write('%g\n' % r_list[j])

        print i, r_list, t_list
        outfile_t.flush()
        #[outfile_r.flush() for outfile_r in outfile_r_list]


    outfile_t.close()
    #[outfile_r.close() for outfile_r in outfile_r_list]



def singlerun(T_list, D_factor, N_B, N_X):

    # 100 nM = 100e-9 * N_A * 100 / m^3 = 6.02e19
    # V = 1 / 6.02e19 = 1.66e-20 m^3
    # L = 2.55e-7 m

    # 1 uM = 6.02e20 / m^3
    # V = 1.66e-21 m^3
    # L = 1.18e-7

    # ### Define constants
    DX_factor = 1

    V = 1e-18 # m^3
    L = V ** (1.0/3.0) 

    matrix_size = min(max(3, int((9 * (N_X+N_B)) ** (1.0/3.0))), 60)
    print 'matrix_size=', matrix_size

    radius = 2.5e-9
    sigma = radius * 2
    r0 = sigma
    D = 1e-12 * D_factor
    D_tot = D * 2

    tau = sigma**2 / D_tot

    #kf = 1000 * sigma * D_tot
    # 1e9 [1 / (M s)] -> 1e9 / 1000 / N_A [m^3 / s]
    kf = 0.092e-18

    DX = D * DX_factor

    # ### Define model
    
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

    # ### Set up simulator

    w = gfrdbase.create_world(L, matrix_size)
    nrw = gfrdbase.create_network_rules_wrapper(m)
    s = EGFRDSimulator(w, myrandom.rng, nrw)

    # ### Place particles 

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


    s.reset()
    # Removed introduction of network rules here.
    # s.set_model(m) # IF presence of network rules results in data
                     # that is not cleaned in reset, this might give 
                     # problems. Originally, s.set_model(m) would 
                     # have introduced these rules here.

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
    """
    while 1:
        pp = s.get_particles_within_radius(A_pos, float(A['radius']))
        if not pp:
            break
        for p in pp:
            s.remove_particle(p)
        s.throw_in_particles(X, len(pp), box1)

    s.place_particle(A, A_pos)

    # Idem for a particle B:
    while 1:
        pp = s.get_particles_within_radius(B_pos, float(B['radius']))
        if not pp:
            break
        for p in pp:
            s.remove_particle(p)
        s.throw_in_particles(X, len(pp), box1)    

    s.place_particle(B, B_pos)

    # Place rest of B at random positions
    if N_B > 1:
        s.throw_in_particles(B, N_B-1, box1)

    r_list = []
    t_list = []
    t_last = 0

    s.step()

    next_stop = T_list[0] # In this case T_list[0] equals infinity

    i_T = 0

    #TODO (following correct?-:)
    """
    What happens here: 
    - If there was a reaction
        - And there are no C-particles, set t_last to current simulation
          time
        - And there are C particles, calculate dt between previous??? and 
          this reaction. 

    OR:
    
        - If the last C particle has just reacted away, record the time
        - If there are still C particles, record the time relative to this  
          occurence 

    T-list in principle contains a list with scheduled "stops"

    """
    while 1:
        if s.last_reaction: # aka if there was a reaction
            print s.last_reaction
            if len(s.world.get_particle_ids(C.id)) == 0:  #A,B
                print 'set t_last', s.t
                t_last = s.t  # set t_last
            else:    # C
                print 'reaction: ', s.t - t_last
                t_list.append(s.t - t_last)

        next_time = s.get_next_time()
        if next_time > next_stop:
            print 'stop', i_T, next_stop
            s.stop(next_stop)
            if len(s.world.get_particle_ids(C.id)) != 0:  #A,B
                r_list.append(0)
            else:
                r_list.append(s.distance_between_particles(A.id, B.id))

            i_T += 1
            next_stop = T_list[i_T]
        
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
