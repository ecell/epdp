
import math
import numpy
import scipy
import scipy.optimize
import myrandom

import _gfrd

Pi = scipy.pi
Pi2 = scipy.pi * 2.0
PiSqrt = math.sqrt(scipy.pi)

N_A = 6.0221367e23
INF = numpy.inf

ZEROPOS = numpy.array([0., 0., 0.])
NOWHERE = numpy.array((INF, INF, INF))

SAFETY = 1.0 + 3e-2     # Lengths of shell in construction are divided by this safety factor 
                        # to make sure that they do not overlap due to numerical rounding errors

# Tolerance used for float comparison functions. Oversimplifying: two floats a 
# and b are considered to be equal if abs(a - b) < TOLERANCE * abs(a).
TOLERANCE = 1e-7
TIME_TOLERANCE = 1e-10

# Scheduler digits: the simulator will round the next-event time steps to this
# precision before putting them into the scheduler. This was originally introduced
# to overcome nasty divergence effects in reloaded simulations caused by limited
# Python float precision. Put this at a high value if you want to sample accurately
# very small time steps by purpose.
# Note that a low number of digits may cause cascades of subsequent dt=0 updates
# for example in pair formation which ultimately will lead to abortion of the simulation.
SCHEDULER_DIGITS = 2*int(-1.0*numpy.log10(TIME_TOLERANCE))

# Multiplication factor used for seperating 2 particles or a particle and a 
# surface after unbinding.
MINIMAL_SEPARATION_FACTOR = 1.0 + TOLERANCE


MULTI_SHELL_FACTOR = math.sqrt(3)
                     # This factor multiplied with the particle radius decides when to add
                     # NonInteractionSingles to a Multi and also defines the Multi shell size.
                     # IMPORTANT NOTE: MULTI_SHELL_FACTOR should be AT LEAST sqrt(2) !
                     # This stems from the fact that there is vacant space in the cylinder

SINGLE_SHELL_FACTOR = 3.5 #2.0*MULTI_SHELL_FACTOR
                      # This is the threshold for when the algorithm switches from forming
                      # NonInteractionSingles to forming a Pair or Interaction. It also defines
                      # the radius in which the NonInteractionSingle will burst intruding domains.
                      # IMPORTANT NOTE: SINGLE_SHELL_FACTOR should be AT_LEAST 2 * MULTI_SHELL_FACTOR !
                      # Otherwise we risk the construction of situations in which two neighbouring
                      # particles neither can burst nor form a multi nor form a minimal single!

CYLINDER_R_FACTOR = 1.0
                    # This factor rescales all cylindrical domains constructed on cylindrical
                    # surfaces. It should be ALWAYS at least 1.0. Keep it at 1.0 if you do not
                    # know what it does.

if __debug__:
    PRECISION = 7
    FORMAT_DOUBLE = '%.' + str(PRECISION) + 'g'


# Float comparison functions.
def feq(a, b, typical=1, tolerance=TOLERANCE):
    """Return True if a and b are equal, subject to given tolerances.  
    Float comparison.

    Also see numpy.allclose().

    The (relative) tolerance must be positive and << 1.0

    Instead of specifying an absolute tolerance, you can specify a 
    typical value for a or b. The absolute tolerance is then the 
    relative tolerance multipied by this typical value, and will be 
    used when comparing a value to zero. By default, the typical 
    value is 1."""

    return abs(a - b) < tolerance * (typical + min(abs(a), abs(b)))

def all_feq(a, b, typical=1, tolerance=TOLERANCE):
    """ Return True if array-like objects a and b are equal,
        subject to given tolerances.  

        Float comparison.

        This is basically a wrapper around numpy.allclose().
    """
    
    return numpy.allclose(a, b, rtol=TOLERANCE, atol=TOLERANCE*typical)
    

def fgreater(a, b, typical=1, tolerance=TOLERANCE):
    """Return True if a is greater than b, subject to given tolerances.  
    Float comparison."""

    return a - b > tolerance * (typical + min(abs(a), abs(b)))


def fless(a, b, typical=1, tolerance=TOLERANCE):
    """Return True if a is less than b, subject to given tolerances.  
    Float comparison."""

    return b - a > tolerance * (typical + min(abs(a), abs(b)))


def fgeq(a, b, typical=1, tolerance=TOLERANCE):
    """Return True if a is greater or equal than b, subject to given 
    tolerances. Float comparison."""

    diff = a - b
    barrier = tolerance * (typical + min(abs(a), abs(b)))
    # Try both 'greater than' and equality.
    return diff > barrier or abs(diff) < barrier


def fleq(a, b, typical=1, tolerance=TOLERANCE):
    """Return True if a is less than or equal than b, subject to given 
    tolerances. Float comparison."""

    diff = b - a
    barrier = tolerance * (typical + min(abs(a), abs(b)))
    # Try both 'less than' and equality.
    return diff > barrier or abs(diff) < barrier

def per_M_to_m3(rate):
    """Convert a reaction rate from units 'per molar per second' to 
    units 'meters^3 per second'.

    """
    return rate / (1000 * N_A)

def per_microM_to_m3(rate):
    """Convert a reaction rate from units 'per micromolar per second' to 
    units 'meters^3 per second'.

    """
    return per_M_to_m3(rate * 1e6)

def M_to_per_m3(molar):
    """Convert a concentration from units 'molar' to units 'per 
    meters^3'.

    """
    return molar * (1000 * N_A)

def microM_to_per_m3(micromolar):
    """Convert a concentration from units 'micromolar' to units 'per 
    meters^3'.

    """
    return M_to_per_m3(micromolar / 1e6) 

def mean_arrival_time(r, D):
    return (r * r) / (6.0 * D)

def uniq(l):
    nset = {}
    map(nset.__setitem__, l, [])
    return nset.keys()

cyclic_transpose = _gfrd.cyclic_transpose

def distance_sq_array_simple(position1, positions, fsize = None):
    return numpy.square(positions - position1).sum(1)

def distance_array_simple(position1, positions, fsize = None):
    return numpy.sqrt(distance_sq_array_simple(position1, positions))

distance = _gfrd.distance

distance_cyclic = _gfrd.distance_cyclic

def distance_sq_array_cyclic(position1, positions, fsize):
    diff = numpy.abs(positions - position1)
    diff -= numpy.greater(diff, fsize * 0.5) * fsize # transpose
    return numpy.square(diff).sum(1)

def distance_array_cyclic(position1, positions, fsize = 0):
    return numpy.sqrt(distance_sq_array_cyclic(position1, positions, fsize))

def cartesian_to_spherical(c):
    # x, y, z = c
    r = length(c)
    theta = math.acos(c[2] / r)
    phi = math.atan2(c[1], c[0])
    if phi < 0.0:  # atan2 returns [- PI, PI]
        phi += 2.0 * Pi
    return numpy.array([r, theta, phi])

def spherical_to_cartesian(s):
    #FIXME: it's possible that the below is a source of some bias. (why??)
    r, theta, phi = s
    # theta is the angle with zenith axis z, phi is the angle with the azimuthal axis x in the direction of the y axis
    # -Pi < theta < Pi
    # -Pi/2 < phi < Pi/2 although since phi is usually a free parameter and drawn from a uniform distribution
    #                    it can also be 0 < phi < 2Pi or 0 < phi < Pi
    sintheta = math.sin(theta)
    return numpy.array([r * math.cos(phi) * sintheta,
                        r * math.sin(phi) * sintheta,
                        r * math.cos(theta)])

def random_unit_vector_s():
    s = numpy.array([1.0, myrandom.uniform(0, Pi), myrandom.uniform(0, Pi2)])
    return s

def random_unit_vector():
    v = [myrandom.uniform(-1,1), myrandom.uniform(-1,1), myrandom.uniform(-1,1)]
    return _gfrd.normalize(v, 1)   

def random_vector(r):
    # v, s = [0.0, 0.0, 0.0], 2.0
    # while s > 1.0:
    #     v[0] = myrandom.uniform (-1, 1)
    #     v[1] = myrandom.uniform (-1, 1)
    #     s = v[0] * v[0] + v[1] * v[1]

    # v[2] = -1 + 2 * s
    # a = 2 * numpy.sqrt(1 - s)
    # v[0] *= a
    # v[1] *= a

    cos_theta = myrandom.uniform(-1, 1)
    sin_theta = numpy.sqrt(1 - cos_theta * cos_theta)
    phi = myrandom.uniform(0, Pi2)
    sin_phi, cos_phi = math.sin(phi), math.cos(phi) # sincos
    v = [sin_theta * cos_phi, sin_theta * sin_phi, cos_theta]
    return _gfrd.normalize(v, r)

def random_vector2D(r):
    # Return a random 2D cartesian vector of length r.
    phi = myrandom.uniform(0, Pi2)
    v = [math.cos(phi), math.sin(phi)] # sincos

    # Todo. return _gfrd.normalize(v, r)
    v = numpy.array(v)
    norm = numpy.linalg.norm(v)
    return v * (r / norm)

def random_sign():
    # NOTE Do not use the numpy.sign function,
    # because it returns 0 for 0 input!
    r = myrandom.uniform(-1,1)
    if r<0:
      return -1
    else:
      return 1
    
def length(a):
    return _gfrd.length(a)

def normalize(a, l=1):
    return _gfrd.normalize(a, l)

def vector_angle(a, b):
    # returns the angle between vector a en b
    # the range of the angle is (0 - Pi)
    cosangle = numpy.dot(a, b) / (length(a) * length(b))
    return math.acos(cosangle)

def vector_angle_against_z_axis(b):
    cosangle = b[2] / length(b)
    return math.acos(cosangle)

def crossproduct(a, b):
    M = numpy.array([[   0.0, - a[2],   a[1]],
                     [  a[2],    0.0, - a[0]],
                     [- a[1],   a[0],    0.0]])
    return numpy.dot(M, b)

def crossproduct_against_z_axis(a):
    # this is crossproduct (a, z-axis)
    return numpy.array([- a[1], a[0], 0.0])

def rotate_vector(v, r, alpha):
    """v: vector to rotate
    r: normalized rotation axis
    alpha: rotation angle in radian"""

    cosalpha = math.cos(alpha)
    sinalpha = math.sin(alpha)
    cosalphac = 1.0 - cosalpha

    M = numpy.array([[cosalpha + cosalphac * r[0] * r[0],
                      cosalphac * r[0] * r[1] - r[2] * sinalpha,
                      cosalphac * r[0] * r[2] + r[1] * sinalpha],
                     [cosalphac * r[0] * r[1] + r[2] * sinalpha,
                      cosalpha + cosalphac * r[1] * r[1],
                      cosalphac * r[1] * r[2] - r[0] * sinalpha],
                     [cosalphac * r[0] * r[2] - r[1] * sinalpha,
                      cosalphac * r[1] * r[2] + r[0] * sinalpha,
                      cosalpha + cosalphac * r[2] * r[2]]])

    return numpy.dot(M,v)

def calculate_pair_CoM(pos1, pos2, D1, D2, world_size):
    return _gfrd.calculate_pair_CoM(pos1, pos2, D1, D2, world_size);

apply_boundary = _gfrd.apply_boundary

def findroot(f, a, b):
    return scipy.optimize.brentq(f, a, b)

def permutate(seq):
    # permutate a sequence and return a list of the permutations
    if not seq:
        return [seq] # is an empty sequence
    else:
        temp = []

        for k in range(len(seq)):
            part = seq[:k] + seq[k+1:]
            for m in permutate(part):
                temp.append(seq[k:k+1] + m)
        return temp

def k_D(Dtot, sigma):
    """Calculate the 'pseudo-'reaction rate (kD) caused by diffusion.
    
    kD is equal to 1 divided by the time it takes for two particles to 
    meet each other by diffusion. It is needed when converting from 
    an intrinsic reaction rate to an overall reaction rates or vice 
    versa.

    Example:
        - A + B -> C.

    Arguments:
        - Dtot:
            the diffusion constant of particle A plus the diffusion 
            constant of particle B. Units: meters^2/second.
        - sigma
            the radius of particle A plus the radius of particle B. 
            Units: meters.

    This function is only available for reaction rules in 3D. No 
    analytical expression for kD in 1D or 2D is currently known. 

    """
    return 4.0 * numpy.pi * Dtot * sigma

def k_a(kon, kD):
    """Convert an overall reaction rate (kon) for a binding/annihilation 
    reaction rule to an intrinsic reaction rate (ka).

    Example:
        - A + B -> C
            binding reaction rule
        - A + B -> 0
            annihilation reaction rule

    Arguments:
        - kon
            the overall reaction rate for the reaction rule. Units: 
            meters^3/second.
        - kD
            the 'pseudo-'reaction rate caused by the diffusion of 
            particles A and B. See the function k_D(). Units: 
            meters^3/second.

    This function is only available for reaction rules in 3D. No 
    analytical expression for kD in 1D or 2D is currently known. 

    """
    if kon > kD:
        raise RuntimeError, 'kon > kD.'
    ka = 1. / ((1. / kon) - (1. / kD))
    return ka

def k_d(koff, kon, kD):
    """Convert an overall reaction rate (koff) for an unbinding reaction 
    rule to an intrinsic reaction rate (kd).

    This one is a bit tricky. We consider reaction rules with only 1 
    reactant. In case there is only 1 product also, no conversion in 
    necessary. But when the reaction rule has 2 products, we need to 
    take the reverse reaction rule into account and do the proper 
    conversion.

    Example:
        - C -> A + B
            unbinding reaction rule
        - A + B -> C
            reverse reaction rule

    Arguments:
        - koff
            the overall reaction rate for the unbinding reaction rule.  
            Units: meters^3/second.
        - kon
            the overall reaction rate for the reverse reaction rule. 
            Units: meters^3/second.
        - kD
            the 'pseudo-'reaction rate caused by the diffusion of 
            particles A and B. See the function k_D(). Units: 
            meters^3/second.

    This function is only available for reaction rules in 3D. No 
    analytical expression for kD in 1D or 2D is currently known. 

    """
    ka = k_a(kon, kD)
    kd = k_d_using_ka(koff, ka, kD)
    return kd

def k_d_using_ka(koff, ka, kD):
    """Convert an overall reaction rate (koff) for an unbinding reaction 
    rule to an intrinsic reaction rate (kd).

    Similar to the function k_d(), but expects an intrinsic rate (ka) 
    instead of an overall rate (kon) for the reversed reaction rule as 
    the second argument.

    This function is only available for reaction rules in 3D. No 
    analytical expression for kD in 1D or 2D is currently known. 

    """
    kd =  koff * (1 + float(ka) / kD)
    return kd

def k_on(ka, kD):
    """Convert an intrinsic reaction rate (ka) for a binding/annihilation 
    reaction rule to an overall reaction rate (kon).

    The inverse of the function k_a().
    
    Rarely needed.

    This function is only available for reaction rules in 3D. No 
    analytical expression for kD in 1D or 2D is currently known. 

    """
    kon = 1. / ((1. / kD) + (1. / ka))  # m^3/s
    return kon

def k_off(kd, kon, kD):
    """Convert an intrinsic reaction rate (kd) for an unbinding reaction 
    rule to an overall reaction rate (koff).

    The inverse of the function k_d().

    Rarely needed.

    This function is only available for reaction rules in 3D. No 
    analytical expression for kD in 1D or 2D is currently known. 

    """
    ka = k_a(kon, kD) 
    koff = k_off_using_ka(kd, ka, kD)
    return koff

def k_off_using_ka(kd, ka, kD):
    """Convert an intrinsic reaction rate (kd) for an unbinding reaction 
    rule to an overall reaction rate (koff).

    Similar to the function k_off(), but expects an intrinsic rate 
    (ka) instead of an overall rate (kon) as the second argument.

    Rarely needed.

    This function is only available for reaction rules in 3D. No 
    analytical expression for kD in 1D or 2D is currently known. 

    """
    koff = 1. / (float(ka) / (kd * kD) + (1. / kd))
    return koff

def C2N(c, V):
    """Calculate the number of particles in a volume 'V' (dm^3) 
    with a concentration 'c' (mol/dm^3).

    """
    return c * V * N_A  # round() here?

