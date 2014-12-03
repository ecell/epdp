# Simulation of the reactions: A + B -(ka)-> C C -(kd)-> A + B
# ### Importing necessary files ###

# import egfrd libraries
from egfrd import * 
import model 
import gfrdbase 
import _gfrd 

# Extra functions
def check_existance( species_id ):
    pids = w.get_particle_ids(species_id)
    if(len(pids) > 0): return True
    else: return False 

def give_pid( species_id ):
    pids = w.get_particle_ids(species_id)
    for pid in pids:
        return pid

# Define the simulation constants.
def constants():
    global L, sigma, D3, ka, kd, Kd, NA, NB, NC_max, rnd_seed, Nstep, tsim
    ka = 1e-19 # Association constant
    kd = 1e3 # Complex dissociation constant
    Kd = kd / ka
    D3 = 1e-12 # Diffusion constant (3D)
    sigma = 1e-9 # Diameter particle
    L = 1e-7 # Lengths of simulation
    NA = 50 # Starting number of A particles
    NB = 50 # Starting number of B particles
    NC_max = 50 #Minimal number of C particles
    rnd_seed = 42 # seed for randomizer
    Nstep = int(1e8) # Number of steps simulation will last
    tsim = 100 # Simulated time
        
    
# Define the model (m) with three particle species (A,B,C) and two reaction rules (r1,r2).
def setup_model():
    global m, A, B, C, r1, r2
    
    # Creates an instance of the ParticleModel class
    m = model.ParticleModel(L)
    # Define species:
    A = model.Species('A', D3, sigma ) # Species([name], [D], [radius=0,]
    B = model.Species('B', D3, sigma )
    C = model.Species('C', D3*D3, sigma )
    # Couple these species to the previously defined particle model
    m.add_species_type(A) # add_species_type(species)
    m.add_species_type(B)
    m.add_species_type(C)
    # Adding reaction rules for B(DNA) + PROMOTOR <-> C(DNA)
    r1 = model.create_binding_reaction_rule(A, B, C, ka)
    m.network_rules.add_reaction_rule(r1)
    
    r2 = model.create_unbinding_reaction_rule(C, A, B, kd)
    m.network_rules.add_reaction_rule(r2)

# SetUp the simulator (s), given the world (w) and the model (m).
def setup_simulator():
    global w, s

    # ### Setting up the simulator: ###
    # Create the world with 3x3x3 cells for the matrixsize
    w = gfrdbase.create_world(m, 3)
    myrandom.seed( rnd_seed )
    # Create a simulator for the world:
    s = EGFRDSimulator(w, myrandom.rng) # EGFRDSimulator(self, world,
                                        # rng, network_rules=None)
    # Place reactant molecules of each species in our world:
    #throw_in_particles(w, A, NA )
    #throw_in_particles(w, B, NB )
    throw_in_particles(w, C, NC_max)

# Defines the number of C-particles given NA, NB, V and Kd
def N_c(Na, Nb):
    part = Na + Nb + Kd * L**3
    NC_min = part - (part**2 - 4 * Na * NB)**0.5
    return NC_min/(2 * L**3)


# ### Set everything up: ###
constants()
setup_model() 
setup_simulator() 
print '* Initialization complete, starting simulation.'

# ### Let's roll: ###
print "NC_min_analytic = %s" %(N_c(NA, NB)) 
t_bind = 0.0 
t_decay = 0.0 
oldt_decay = 0.0 
Nbinds = 0 
step_cntr = 0 
lastC_pid = give_pid( A ) 
CatLastReaction = False; 
tCbound = 0.0 

while(s.get_next_time() < tsim):
    s.step()
    if step_cntr%10000 == 1000:
        print "step#: %s mln, t = %s, #rebinds = %s, #A = %s, #B = %s, #C= %s, Bound fraction = %s, [C]= %s" %(step_cntr/1e6, s.t, Nbinds, len(w.get_particle_ids(A)), len(w.get_particle_ids(B)), len(w.get_particle_ids(C)), tCbound/s.t, len(w.get_particle_ids(C))/L**3)
  
    if s.last_reaction:
        #Check is last_reaction was the formation of a NEW C particle.
        if check_existance( C ) and give_pid( C ) != lastC_pid:
            lastC_pid = give_pid( C )
            t_bind = s.t
            Nbinds += 1
        #Check if last_reaction was the dissociation of a C particle.
        if CatLastReaction and not check_existance( C ):
            oldt_decay = t_decay
            t_decay = s.t
            tCbound += t_decay - t_bind
            
        CatLastReaction = check_existance( C )
    
    step_cntr += 1

# Synchronize all domains to tsim.
s.stop(tsim)

