#!/usr/bin/python

# This script tries to test the GF in the whole parameter regime.






# Modules
# ===============================
# Standard stuff
import sys
import os
import datetime

from egfrd import *
#from visualization import vtklogger

import model
import gfrdbase
import _gfrd

# GF itself
# ===
import _greens_functions



# Some typical "crash" Constants
# ===============================
d = 0.931552; r = 6.10119e-08; dt = 0.00114661; D = 2e-12; sigma = 2e-09; a = 6.10119e-08; 
kf = 0; r0 = 4.95059e-09; h = 0; rnd = 0.24;

# Let's get some data 
# ===============================

mygf = _greens_functions.GreensFunction2DRadAbs(D, kf, r0, sigma, r)
try:
    mygf.drawTheta(rnd, r, dt)
    print "succes!"   
except:
    print "fail!"



# Let's try and NOT make it crash
# ===============================
D = 2e-12; 
sigma = 2e-09; # 
a = sigma*5; ##
r = sigma*3; ## 
r0 = sigma*2; ## 
dt = 0.03; ##
kf = 0; ##
h = 0; ##

rnd = 0.1; ##

# Let's get some data 
# ===============================

mygf = _greens_functions.GreensFunction2DRadAbs(D, kf, r0, sigma, r)
mygf.drawTheta(rnd, r, dt)

try:
    mygf.drawTheta(rnd, r, dt)
    print "succes!"   
except:
    print "fail!"























