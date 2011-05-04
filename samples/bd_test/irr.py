#!/usr/bin/env python

from bd import *
import sys
import gfrdbase
import model

def run(outfilename, T, N):

    outfile = open(outfilename, 'w')

    for i in range(N):
        d, t = singlerun(T)
        outfile.write('%g\n' % d)

        print d, t
        #assert d == 0 or t == T

    outfile.close()


def singlerun(T):

    sigma = 5e-9
    r0 = sigma
    D = 1e-12
    kf = 10 * sigma * D

    m = model.ParticleModel(world_size=1e-5)

    A = model.Species('A', D/2, sigma/2)
    B = model.Species('B', D/2, sigma/2)
    C = model.Species('C', D/2, sigma/2)
    # Couple these species to the previously defined particle model
    m.add_species_type(A) # "m.add_species_type(species=A)" not
                          # allowed due boost package
    m.add_species_type(B) 
    m.add_species_type(C) 
    
    r1 = model.create_binding_reaction_rule(A, B, C, kf)

    m.network_rules.add_reaction_rule(r1)

    s.set_model(m)
    
    particleA = s.place_particle(A, [0,0,0])
    particleB = s.place_particle(B, [(float(A['radius']) + float(B['radius']))+1e-23,0,0])

    end_time = T

    # Create world size m, with matrix_size sized
    # neighbour matrix
#    w = gfrdbase.create_world(m=1e-3, matrix_size=3)
    w = gfrdbase.create_world(1e-3, 3)
    s = BDSimulator(w)


    #s.initialize()

    while 1:
        next_time = s.t + s.dt
        if next_time > end_time:
            break
        s.step()

        if s.core.last_reaction:
            print 'reaction'
            return 0.0, s.t

    distance = s.distance_between_particles(particleA[1], particleB[1])

    return distance, s.t


def first(x):
    x = iter(x)
    try:
        return x.next()
    except StopIteration, e:
        return None
    

if __name__ == '__main__':
    run(sys.argv[1], float(sys.argv[2]), int(sys.argv[3]))
