#!/usr/env python

import math
import numpy

from single import *
from pair import *
from multi import *

__all__ = [ 'Histogram3D', 'DomainsHistogramCollection' ]


class Histogram3D(object):

    def __init__(self, nbins, dimensions):

        self.nbins      = nbins
        self.dimensions = dimensions
        self.histogram  = numpy.zeros(nbins)
        self.total_cnt  = 0

        self.binsizes   = self.calculate_binsizes()

        self.normalized = False


    def calculate_binsizes(self):

        binsizes = numpy.zeros(3)

        for d in [0, 1, 2]:

            binsizes[d] = 1.0 * self.dimensions[d] / self.nbins[d]

        return binsizes


    def bin(self, data):

        index = [-1, -1, -1]

        for d in [0, 1, 2]:

            b = int(data[d] / self.binsizes[d])

            if b > 0 and b < self.nbins[d]-1:

                    index[d] = b

        ilist = [index[0], index[1], index[2]]

        if all([i>=0 and i<=self.nbins for i in ilist]):
            
            [ix, iy, iz] = ilist
            self.histogram[ix][iy][iz] = self.histogram[ix][iy][iz] + 1

            self.total_cnt = self.total_cnt + 1


    def normalize(self):
        # TODO
        pass   
    

class DomainsHistogramCollection(object):

    def __init__(self, world, nbins):

        self.world = world
        self.ws = world.world_size

        self.nbins = nbins

        self.SinglesHistogram   = Histogram3D([self.nbins, self.nbins, self.nbins], [self.ws, self.ws, self.ws])
        self.PairsHistogram     = Histogram3D([self.nbins, self.nbins, self.nbins], [self.ws, self.ws, self.ws])
        self.MultisHistogram    = Histogram3D([self.nbins, self.nbins, self.nbins], [self.ws, self.ws, self.ws])
        self.NonMultisHistogram = Histogram3D([self.nbins, self.nbins, self.nbins], [self.ws, self.ws, self.ws])


    def bin_domain(self, domain):

        if isinstance(domain, Single):

            pos = domain.pid_particle_pair[1].position

            self.SinglesHistogram.bin(pos)
            self.NonMultisHistogram.bin(pos)

        elif isinstance(domain, Pair):

            pos1 = domain.pid_particle_pair1[1].position
            pos2 = domain.pid_particle_pair2[1].position

            for pos in [pos1, pos2]:
                self.PairsHistogram.bin(pos)
                self.NonMultisHistogram.bin(pos)

        elif isinstance(domain, Multi):

            for pid_particle_pair in domain.particle_container:

                pos = pid_particle_pair[1].position            
                self.MultisHistogram.bin(pos)


    def normalize(self):
        # TODO
        pass

