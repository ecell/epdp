#!/usr/bin/env python

from egfrd import *

from logger import *
import sys

import math

#V = 1e-15 #liter
#L = 1e-6
V = 3.33e-17
L = math.pow( V, 1.0 / 3.0 )

s = EGFRDSimulator()
s.setCellSize( L )

box1 = CuboidalSurface( [0,0,0],[L,L,L] )
# not supported yet
#s.addSurface( box1 )

K = Species( 'K', 2e-12, 5e-9 )
s.addSpecies( K )
KK = Species( 'KK', 2e-12, 5e-9 )
s.addSpecies( KK )
P = Species( 'P', 2e-12, 5e-9 )
s.addSpecies( P )
Kp = Species( 'Kp', 2e-12, 5e-9 )
s.addSpecies( Kp )
Kpp = Species( 'Kpp', 2e-12, 5e-9 )
s.addSpecies( Kpp )
K_KK = Species( 'K_KK', 2e-12, 5e-9 )
s.addSpecies( K_KK )
Kp_KK = Species( 'Kp_KK', 2e-12, 5e-9 )
s.addSpecies( Kp_KK )
Kpp_KK = Species( 'Kpp_KK', 2e-12, 5e-9 )
s.addSpecies( Kpp_KK )
Kpp_P = Species( 'Kpp_P', 2e-12, 5e-9 )
s.addSpecies( Kpp_P )
Kp_P = Species( 'Kp_P', 2e-12, 5e-9 )
s.addSpecies( Kp_P )

#  1 2   K + KK   <-> K_KK
#  3     K_KK       -> Kp + KK
#  4 5   Kp + KK  <-> Kp_KK
#  6     Kp_KK      -> Kpp + KK 
#  7 8   Kpp + P <-> Kpp_P
#  9     Kpp_P     -> Kp + P
# 10 11  Kp + P  <-> Kp_P
# 12     Kp_P      -> K + P
# 12     Kp_P      -> K + P
# 13     Kpp     -> Kp
# 14     Kp      -> K


Dtot = 4e-12
sigma = 5e-9
Dpisigma4 = 4 * numpy.pi * Dtot * sigma

def k_i( k ):
    kk = k * 1e3 / N_A
    return - kk * Dpisigma4 / (kk - Dpisigma4)

print k_i( 0.02e-9)


r1 = BindingReactionType( K, KK, K_KK, k_i(0.02e-9) )
s.addReactionType( r1 )
r2 = UnbindingReactionType( K_KK, K, KK, 1 )
s.addReactionType( r2 )
r3 = UnbindingReactionType( K_KK, Kp, KK, 1e-2 )
s.addReactionType( r3 )

r4 = BindingReactionType( Kp, KK, Kp_KK, k_i( 0.032e-9 ) )
s.addReactionType( r4 )
r5 = UnbindingReactionType( Kp_KK, Kp, KK, 1 )
s.addReactionType( r5 )
r6 = UnbindingReactionType( Kp_KK, Kpp, KK, 15 )
s.addReactionType( r6 )

r7 = BindingReactionType( Kpp, P, Kpp_P, k_i( 0.02e-9 ) )
s.addReactionType( r7 )
r8 = UnbindingReactionType( Kpp_P, Kpp, P, 1 )
s.addReactionType( r8 )
r9 = UnbindingReactionType( Kpp_P, Kp, P, 0.01 )
s.addReactionType( r9 )

r10 = BindingReactionType( Kp, P, Kp_P, k_i( 0.032e-9 ))
s.addReactionType( r10 )
r11 = UnbindingReactionType( Kp_P, Kp, P, 1 )
s.addReactionType( r11 )
r12 = UnbindingReactionType( Kp_P, K, P, 15 )
s.addReactionType( r12 )

#r13 = UnimolecularReactionType( Kpp, Kp, 1e-1 )
#s.addReactionType( r13 )
#r14 = UnimolecularReactionType( Kp, K, 1e-1 )
#s.addReactionType( r14 )

def C2N( c ):
    return round( c * V * N_A )

s.throwInParticles( K, C2N( 500e-9 ), box1 )

s.throwInParticles( KK, C2N( 50e-9 ), box1 )
s.throwInParticles( P, C2N( 50e-9 ), box1 )


l = Logger( s, 'mapk' )
#l.setParticleOutput( ('Ea','X','EaX','Xp','Xpp','EaI') )
l.setInterval( 1e-0 )
l.log()


while s.t < 100:
    s.step()
    s.dumpPopulation()
    l.log()
    

