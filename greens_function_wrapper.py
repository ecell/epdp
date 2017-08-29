import numpy
import myrandom
from _greens_functions import PairEventKind
from constants import EventType
from utils import *

import logging
log = logging.getLogger('ecell')


def draw_time_wrapper(gf):
    rnd = myrandom.uniform()

    if __debug__:
        log.debug('        * drawTime: ' + gf.__class__.__name__)
    try:
        dt = gf.drawTime(rnd)
    except Exception, e:
        raise Exception('gf.drawTime() failed, '
                        '%s, rnd = %g, %s' %
                        (str(e), rnd, gf.dump()))
    return dt

def draw_event_type_wrapper(gf, dt):
    rnd = myrandom.uniform()

    if __debug__:
        log.debug('        * drawEventType: ' + gf.__class__.__name__)
    try:
        event_type = gf.drawEventType(rnd, dt)
    except Exception, e:
        raise Exception('gf.drawEventType() failed, '
                        '%s, rnd = %g, dt = %g, %s' %
                        (str(e), rnd, dt, gf.dump()))
    return event_type

def draw_r_wrapper(gf, dt, a, sigma=None):
    rnd = myrandom.uniform()

    if __debug__:
        log.debug('        * drawR: ' + gf.__class__.__name__)
    try:
        r = gf.drawR(rnd, dt)
        while r > a or r <= sigma: # redraw; shouldn't happen often
            if __debug__:
                log.debug('        * drawR: redraw')
            rnd = myrandom.uniform()
            r = gf.drawR(rnd, dt)
    except Exception, e:
        raise Exception('gf.drawR() failed, '
                        '%s, rnd = %g, dt = %g, %s' %
                        (str(e), rnd, dt, gf.dump()))

    return r

def draw_theta_wrapper(gf, r, dt):
    """Draw theta for the inter-particle vector.

    """
    # We're interested in value between [-pi,pi], but this problem is symmetric: 
    # the cumulative pdf is odd. Thus drawTheta() returns a value between 
    # [0,pi], which is then with a 50% probability made negative or kept 
    # positive.
    rnd = myrandom.uniform()

    if __debug__:
        log.debug('        * drawTheta: ' + gf.__class__.__name__)
    try:
        theta = gf.drawTheta(rnd, r, dt)
    except Exception, e:
        print 'gf.drawTheta() failed, %s, rnd = %g, r = %g, dt = %g; %s' %(str(e), rnd, r, dt, gf.dump())
        #ugly hack: When drawTheta fails, return uniformly distributed theta.
        theta = rnd * numpy.pi

    # Heads up. For cylinders theta should be between [-pi, pi]. For 
    # spheres it doesn't matter.
    return myrandom.choice(-1, 1) * theta

