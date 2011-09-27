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




# List of parameters to investigate:
# ======================================

D =      2e-12
sigma =  2e-9
a =      [10*sigma, 100*sigma]
r0 =     [1.001*sigma, 3*sigma, 4.999*sigma]
dt =     [1e-4, 0.03, 0.3, 1, 1e4]
rnd =    [1e-3, 0.05, 0.1, 0.25, (0.5-1e-3)]
kf =     0; # 2e-12 10 1e10]
h =      0;

# Convenient for display later on
# ======================================

codinga = [' - ', ' + ']
codingr0 = [' - ', ' 0 ', ' + ']
codingdt = [' --', '  -', ' 0 ', ' + ', ' ++']
codingrnd = [' --', '  -', ' 0 ', ' + ', ' ++']

resultstable = ['| a | r0| dt|rnd|  -> Result']

# Let's do this the simple way:
# ======================================

if __name__=="__main__":

    print "Starting tests."
    print "========================="

    for i in range(len(a)):
        r = a[i] #!
        for j in range(len(r0)):
            for k in range(len(dt)):
                 for l in range(len(rnd)):


                    print("Parameters: {D = "+str(D)+", sigma = "+
                          str(sigma)+", a = "+str(a[i])+", r0 = "+
                          str(r0[j])+", dt = "+str(dt[k])+", kf = "+
                          str(kf)+", h = "+str(h)+", r = "+str(r)+", "+
                          "rnd = "+str(rnd[l])+"}")
                    
                    mygf = _greens_functions.GreensFunction2DRadAbs(D, kf, r0[j], sigma, r)

                    try:
                        mygf.drawTheta(rnd[l], r, dt[k])
                        print "succes!"   
                        succes = True;
                    except:
                        print "fail!"
                        print sys.exc_info()
                        succes = False;

                    print ""

                    resultstable.append("|"+str(codinga[i])+"|"+
                                    str(codingr0[i])+"|"+
                                    str(codingdt[i])+"|"+
                                    str(codingrnd[i])+
                                    "|  -> "+str(succes)+"")

print "\n\nRESULTS\n========================\n\n"
for a in range(len(resultstable)):
    print str(resultstable[a]);                            
print "\n\n===========================\n\n"


















































