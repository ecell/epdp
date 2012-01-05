from _gfrd import *
from constants import EventType
from _greens_functions import *
from greens_function_wrapper import *

from domain import (
    Domain,
    ProtectiveDomain)
from single import (
    NonInteractionSingle,
    InteractionSingle)
from multi import (
    Multi)
from shells import *
from utils import *

__all__ = [
    'CylindricalSurfacePair',
    'PlanarSurfacePair',
    'SphericalPair',
    'Pair',
    'SimplePair',
    'MixedPair2D3D',
    'MixedPair1D3D',
    ]

if __debug__:
    PRECISION = 3
    FORMAT_DOUBLE = '%.' + str(PRECISION) + 'g'


class Pair(ProtectiveDomain):
    """There are 3 types of pairs:
        * SphericalPair
        * PlanarSurfacePair
        * CylindricalSurfacePair

    """
    # CUTOFF_FACTOR is a threshold to choose between the real and 
    # approximate Green's functions.
    # H = 4.0: ~3e-5, 4.26: ~1e-6, 5.0: ~3e-7, 5.2: ~1e-7,
    # 5.6: ~1e-8, 6.0: ~1e-9
    CUTOFF_FACTOR = 5.6

#    def __init__(self, domain_id, single1, single2, shell_id, r0, 
#                 rrs):
    def __init__(self, domain_id, shell_id, rrs):
        ProtectiveDomain.__init__(self, domain_id, shell_id)

        assert self.testShell       # assume that the testShell exists
        assert isinstance(self.testShell, testPair)
        self.multiplicity = 2

        self.single1 = self.testShell.single1
        self.single2 = self.testShell.single2 
        self.pid_particle_pair1 = self.testShell.pid_particle_pair1
        self.pid_particle_pair2 = self.testShell.pid_particle_pair2

        self.interparticle_rrs = rrs
        self.interparticle_ktot = self.calc_ktot(rrs)
        self.com   = self.testShell.com
        self.iv    = self.testShell.iv
        self.r0    = self.testShell.r0
        self.sigma = self.testShell.sigma

#        self.shell_id = shell_id

    def __del__(self):
        if __debug__:
            log.debug('del %s' % str(self))

    def get_D_tot(self):
        return self.pid_particle_pair1[1].D + \
               self.pid_particle_pair2[1].D
    D_tot = property(get_D_tot)
    D_r   = property(get_D_tot)

    def get_D_R(self):
        return (self.pid_particle_pair1[1].D *
                self.pid_particle_pair2[1].D) / self.D_tot
    D_R = property(get_D_R)

    def initialize(self, t):
        # re-initialize the time parameters of the Domain to reuse it.
        # THIS IS PROBABLY USELESS
        self.last_time = t
        self.dt = 0
        self.event_type = None

    # draws a first event time and the (not fully specified) type of event
    def draw_com_escape_or_iv_event_time_tuple(self, r0):
        """Returns a (event time, event type, reactingsingle=None) tuple.
        
        """
        dt_com = draw_time_wrapper(self.com_greens_function())
        dt_iv = draw_time_wrapper(self.iv_greens_function(r0))  # uses the full solution to draw IV first event time
        if dt_com < dt_iv:
            return dt_com, EventType.COM_ESCAPE, None
        else:
            # Note: we are not calling pair.draw_iv_event_type yet, but 
            # postpone it to the very last minute (when this event is 
            # executed in fire_pair). So IV_EVENT can still be an iv 
            # escape or an iv reaction.
            return dt_iv, EventType.IV_EVENT, None

    def draw_single_reaction_time_tuple(self):
        """Return a (reaction time, event type, reactingsingle)-tuple.

        """
        dt_reaction1, event_type1 = self.single1.draw_reaction_time_tuple()
        dt_reaction2, event_type2 = self.single2.draw_reaction_time_tuple()
        if dt_reaction1 < dt_reaction2:
            return dt_reaction1, event_type1, self.single1
        else:
            return dt_reaction2, event_type2, self.single2

    def determine_next_event(self, r0):
        """Return a (event time, event type, reactingsingle)-tuple.

        """
        return min(self.draw_com_escape_or_iv_event_time_tuple(r0), 
                   self.draw_single_reaction_time_tuple()) 

    def draw_iv_event_type(self, r0):
        gf = self.iv_greens_function(r0)
        event_kind = draw_event_type_wrapper(gf, self.dt)
        if event_kind == PairEventKind.IV_REACTION:
            return EventType.IV_REACTION
        elif event_kind == PairEventKind.IV_ESCAPE:
            return EventType.IV_ESCAPE
        raise NotImplemented()

    def draw_new_positions(self, dt, r0, old_iv, event_type):
        """Calculate new positions of the pair particles using a new 
        center-of-mass, a new inter-particle vector, and an old 
        inter-particle vector.

        """
        new_com = self.draw_new_com(dt, event_type)
        new_iv = self.draw_new_iv(dt, r0, old_iv, event_type)

        unit_z = self.get_shell_direction()
        newpos1, newpos2 = self.do_back_transform(new_com, new_iv,
                                                  self.pid_particle_pair1[1].D,
                                                  self.pid_particle_pair2[1].D,
                                                  self.pid_particle_pair1[1].radius,
                                                  self.pid_particle_pair2[1].radius,
                                                  self.structure,
                                                  unit_z)

        return newpos1, newpos2

    def draw_new_com(self, dt, event_type):
        # draws a new coordinate for the CoM in world coordinates
        if event_type == EventType.COM_ESCAPE:
            r = self.a_R
        else:
            gf = self.com_greens_function()
            r = draw_r_wrapper(gf, dt, self.a_R)

        # Add displacement to old CoM. This assumes (correctly) that 
        # r0=0 for the CoM. Compare this to 1D singles, where r0 is not  
        # necesseraly 0.

        # note that we need to make sure that the com and com_vector are in the structure to prevent
        # the particle leaving the structure due to numerical errors
        return self.com + self.create_com_vector(r)

    def draw_new_iv(self, dt, r0, old_iv, event_type):
        gf = self.choose_pair_greens_function(r0, dt)
        if event_type == EventType.IV_ESCAPE:
            r = self.a_r
        elif event_type == EventType.IV_REACTION:
            r = self.sigma      # maybe this should be zero (product particle is at CoM)
        else:
            r = draw_r_wrapper(gf, dt, self.a_r, self.sigma)

        # note that we need to make sure that the interparticle vector is in the structure to prevent
        # the particle leaving the structure due to numerical errors
        return self.create_interparticle_vector(gf, r, dt, r0, old_iv)

    def check(self):
    # perform internal consistency check 
        pass

    def __str__(self):
        sid1 = self.pid_particle_pair1[1].sid
        sid2 = self.pid_particle_pair2[1].sid
        name1 = self.world.model.get_species_type_by_id(sid1)["name"]
        name2 = self.world.model.get_species_type_by_id(sid2)["name"]
        return 'Pair[%s: %s, %s, (%s, %s)]' % (
            self.domain_id,
            self.pid_particle_pair1[0],
            self.pid_particle_pair2[0],
            name1, name2)

class SimplePair(Pair):

    @classmethod
    def get_min_shell_size(cls, single1, single2, geometrycontainer):
        # 1. Determine min shell size for a putative SimplePair.
        # TODO Note that the calculation is dimension independent. Also
        # 'size' is not the right term here

        assert single1.structure == single2.structure

        pos1 = single1.pid_particle_pair[1].position
        pos2 = single2.pid_particle_pair[1].position
        radius1 = single1.pid_particle_pair[1].radius
        radius2 = single2.pid_particle_pair[1].radius

#        assert (pos1 - single1.shell.shape.position).sum() == 0        # TODO Not sure why this is here
        sigma = radius1 + radius2
        r0 = geometrycontainer.world.distance(pos1, pos2)

        assert r0 >= sigma, \
            'distance_from_sigma (pair gap) between %s and %s = %s < 0' % \
            (single1, single2, FORMAT_DOUBLE % (r0 - sigma))

        # calculate the minimum shell size for the interparticle domain
        D1 = single1.pid_particle_pair[1].D
        D2 = single2.pid_particle_pair[1].D
        D12 = D1 + D2
        dist_from_com1 = r0 * D1 / D12                  # particle distance from CoM
        dist_from_com2 = r0 * D2 / D12
        iv_shell_size1 = dist_from_com1 + radius1       # the shell should surround the particles
        iv_shell_size2 = dist_from_com2 + radius2

        # also fix the minimum shellsize for the CoM domain
        com_shell_size = max(radius1, radius2)

        # calculate total radii including the margin for the burst volume for the particles
        shell_size = max(iv_shell_size1 + com_shell_size + radius1 * (1 - SINGLE_SHELL_FACTOR),
                         iv_shell_size2 + com_shell_size + radius2 * (1 - SINGLE_SHELL_FACTOR))

        return shell_size


    @classmethod
    def get_max_shell_size(cls, shell_center, single1, single2, geometrycontainer, domains):
        # returns the maximum shell size based on geometric constraints
        # It makes sure that surrounding NonInteractionSingles have at least the
        # 'bursting volume' of space.

#        pos1 = single1.pid_particle_pair[1].position
#        pos2 = single2.pid_particle_pair[1].position
#        D1 = single1.pid_particle_pair[1].D
#        D2 = single2.pid_particle_pair[1].D
#        D12 = D1 + D2

        # TODO it makes no sence to calculate the CoM in the world object
#        com  = geometrycontainer.world.calculate_pair_CoM(pos1, pos2, D1, D2)
#        shell_center = geometrycontainer.world.apply_boundary(com)

        # TODO get_closest_obj searches the spherical surrounding although this doesn't make sence for
        # a 2D or 1D pair
#        closest, closest_distance = geometrycontainer.get_closest_obj(shell_center, domains, ignore=[single1.domain_id,
#                                                                      single2.domain_id],
#                                                                      ignores=[single1.structure.id])

        neighbor_domains = geometrycontainer.get_neighbor_domains(shell_center, domains, ignore=[single1.domain_id, 
                                                                  single2.domain_id])
        neighbor_surfaces = geometrycontainer.get_neighbor_surfaces(shell_center, ignores=[single1.structure.id])
        neighbors = neighbor_domains + neighbor_surfaces


        # 3. Check if bursted Singles can still make a minimal shell.
        # TODO this is not correct. If the closest domain is not a NonInteractionSingle
        # we still have to make sure that the (non closest) NonInteractionSingles have
        # the required space

        max_shell_size = geometrycontainer.get_max_shell_size()

        distance_domain_closest = (max_shell_size, None)
        for domain, distance in neighbors:

            if isinstance(domain, NonInteractionSingle):
                closest_particle_distance = geometrycontainer.world.distance(
                        shell_center, domain.pid_particle_pair[1].position)

                distance = min(distance,
                               closest_particle_distance - domain.pid_particle_pair[1].radius * SINGLE_SHELL_FACTOR)

                # options for shell size:
                # a. upto the closest shell assuming that this shell is bigger than the minimum
                # b. The distance to a bursted single including its minimal shell
            distance_domain = (distance, domain)    
            distance_domain_closest = min(distance_domain_closest, distance_domain)

        closest_distance, closest = distance_domain_closest

#                shell_size = min(closest_distance, closest_particle_distance - closest_min_shell)
#                assert shell_size <= closest_distance
#
#            else:
#                assert isinstance(closest, (InteractionSingle, Pair, Multi, Surface, None.__class__))
#                shell_size = closest_distance 

        if __debug__:
            log.debug('Pair closest neighbor: %s, distance: %s' % \
                      (closest, FORMAT_DOUBLE % closest_distance))

        # shortcut
        # if closest_distance <0 then neighbor particle overlaps with center of mass ->
        # neighboring particle is actually closer to the CoM than any of the two member particles
        # this can happen because we also try to make pairs with the next nearest neighbor if the nearest
        # neighbor fails
        if closest_distance <= 0:
            return 0
        else:
            # FIXME Bloody ugly tweak to make sure that in case of a cylinder the 'maximum size' (which is
            # actually the radius) does not stick out of the maximum sphere (which is the real hard max).
            # BLoody UGLY because it now applies to ALL SimplePairs.
            sigma = max(single1.pid_particle_pair[1].radius, single2.pid_particle_pair[1].radius)
            closest_distance = math.sqrt(abs(closest_distance**2 - sigma**2))
            return closest_distance/SAFETY


#    def __init__(self, domain_id, shell_center, single1, single2, shell_id,
#                      r0, shell_size, rrs, structure):
    def __init__(self, domain_id, shell_id, rrs):
        Pair.__init__(self, domain_id, shell_id, rrs)

        assert self.testShell
        assert isinstance(self.testShell, testSimplePair)
        self.structure = self.testShell.structure              # structure on which both particles live

        shell_size = self.get_shell_size()                  # FIXME
        self.a_R, self.a_r = self.determine_radii(self.r0, shell_size)
        # set the radii of the inner domains as a function of the outer protective domain

# this can go
#        self.shell = self.create_new_shell(shell_center, shell_size, domain_id)


#    def getSigma(self):
#        return self.pid_particle_pair1[1].radius + \
#               self.pid_particle_pair2[1].radius
#    sigma = property(getSigma)
#
#    def get_com(self):
#        return self.shell.shape.position
#    com = property(get_com)

#    @classmethod
#    def do_transform(cls, single1, single2, world):
#
#        pos1 = single1.pid_particle_pair[1].position
#        pos2 = single2.pid_particle_pair[1].position
#        D1 = single1.pid_particle_pair[1].D
#        D2 = single2.pid_particle_pair[1].D
#
#        com = world.calculate_pair_CoM(pos1, pos2, D1, D2)
#        com = world.apply_boundary(com)
#        # make sure that the com is in the structure of the particle (assume that this is done correctly in
#        # calculate_pair_CoM)
#
#        pos2t = world.cyclic_transpose(pos2, pos1)
#        iv = pos2t - pos1
#
#        return com, iv

    @ classmethod
    def do_back_transform(cls, com, iv, D1, D2, radius1, radius2, surface, unit_z):
    # here we assume that the com and iv are really in the structure and no adjustments have to be
    # made

        D_tot = D1 + D2
        pos1 = com - iv * (D1 / D_tot)
        pos2 = com + iv * (D2 / D_tot)

        return pos1, pos2

    def get_shell_size(self):
        return self.shell.shape.radius

    def determine_radii(self, r0, shell_size):
        """Determine a_r and a_R from the size of the protective domain.

        Todo. Make dimension (1D/2D/3D) specific someday. Optimization only.

        """
        radius1 = self.pid_particle_pair1[1].radius
        radius2 = self.pid_particle_pair2[1].radius

        D1 = self.pid_particle_pair1[1].D
        D2 = self.pid_particle_pair2[1].D

        LD_MAX = 20 # temporary value
        a_r_max = LD_MAX * (r0 - self.sigma) + self.sigma
        # a_r_max is the maximum size of a_r for a given maximum ratio of l/delta

        # Make sure that D1 != 0 to avoid division by zero in the followings.
        if D1 == 0:
            D1, D2 = D2, D1     # shouldn't we also here swap the particle radii if we're swapping diffusion contants?

        shell_size /= SAFETY    #??

        D_tot = D1 + D2
        D_geom = math.sqrt(D1 * D2)

        assert r0 >= self.sigma, \
            '%s;  r0 %g < sigma %g' % (self, r0, self.sigma)

        # equalize expected mean t_r and t_R.
        if ((D_geom - D2) * r0) / D_tot + shell_size +\
                math.sqrt(D2 / D1) * (radius1 - shell_size) - radius2 >= 0:
            Da = D1
            Db = D2
            radiusa = radius1
            radiusb = radius2
        else:
            Da = D2
            Db = D1
            radiusa = radius2
            radiusb = radius1


        #aR
        a_R = (D_geom * (Db * (shell_size - radiusa) + \
                         Da * (shell_size - r0 - radiusa))) /\
              (Da * Da + Da * Db + D_geom * D_tot)

        #ar
        a_r = (D_geom * r0 + D_tot * (shell_size - radiusa)) / (Da + D_geom)


        # Now if the planned domainsize for r is too large for proper convergence of the Green's
        # functions, make it the maximum allowed size
        if (a_r > a_r_max):
            a_r = a_r_max
            a_R = shell_size - radiusa - a_r * Da / D_tot
            if __debug__:
                log.info('domainsize changed for convergence.')

        assert a_R + a_r * Da / D_tot + radiusa >= \
               a_R + a_r * Db / D_tot + radiusb

        assert abs(a_R + a_r * Da / D_tot + radiusa - shell_size) \
            < 1e-12 * shell_size                        # here the shell_size is the relevant scale


        if __debug__:
          log.debug('a %g, r %g, R %g r0 %g' % 
                 (shell_size, a_r, a_R, r0))
        if __debug__:
            tr = ((a_r - r0)**2) / (6 * self.D_r)       # the expected escape time of the iv
            if self.D_R == 0:
                tR = numpy.inf 
            else:
                tR = (a_R**2) / (6 * self.D_R)          # the expected escape time of the CoM
            log.debug('tr %g, tR %g' % (tr, tR))


        assert a_r > 0
        assert a_r > r0, '%g %g' % (a_r, r0)
        assert a_R > 0 or (a_R == 0 and (D1 == 0 or D2 == 0))

        return a_R, a_r

    def __str__(self):
        pass

#class SphericalPair(SimplePair):
class SphericalPair(SimplePair, Others, hasSphericalShell):
    """2 Particles inside a (spherical) shell not on any surface.

    """
#    def __init__(self, domain_id, shell_center, single1, single2, shell_id,
#                 r0, shell_size, rrs, structure):
    def __init__(self, domain_id, shell_id, testShell, rrs):
#        single1 = testShell.single1
#        single2 = testShell.single2
#        structure = testShell.structure
#        shell_center = testShell.center
#        shell_size = testShell.radius  # TODO this should be removed
#        r0 = testShell.r0        # TODO this should be removed

#        SimplePair.__init__(self, domain_id, shell_center, single1, single2, shell_id,
#                      r0, shell_size, rrs, structure)
        hasSphericalShell.__init__(self, testShell, domain_id)
        SimplePair.__init__(self, domain_id, shell_id, rrs) # Always initialize after hasSphericalShell


    def com_greens_function(self):
        # Green's function for centre of mass inside absorbing sphere.
        return GreensFunction3DAbsSym(self.D_R, self.a_R)

    def iv_greens_function(self, r0):
        # Green's function for interparticle vector inside absorbing 
        # sphere.  This exact solution is used for drawing times.
        return GreensFunction3DRadAbs(self.D_r, self.interparticle_ktot, r0,
                                              self.sigma, self.a_r)

# this can go
#    def create_new_shell(self, position, radius, domain_id):
#        return SphericalShell(domain_id, Sphere(position, radius))

# this is a part of hasShell (see testShell.get_orientation_vector())
    def get_shell_direction(self):
        return random_vector(1.0)

    # selects between the full solution or an approximation where one of
    # the boundaries is ignored
    def choose_pair_greens_function(self, r0, t):
        distance_from_sigma = r0 - self.sigma
        distance_from_shell = self.a_r - r0

        threshold_distance = Pair.CUTOFF_FACTOR * \
            math.sqrt(6.0 * self.D_r * t)

        # if sigma reachable
        if distance_from_sigma < threshold_distance:
        
            # if shell reachable
            if distance_from_shell < threshold_distance:
                # near both a and sigma;
                # use GreensFunction3DRadAbs
                if __debug__:
                    log.debug('GF: normal')
                return self.iv_greens_function(r0)
            else:
                # near sigma; use GreensFunction3DRadInf
                if __debug__:
                    log.debug('GF: only sigma')
                return GreensFunction3DRadInf(self.D_r, self.interparticle_ktot, r0,
                                               self.sigma)

        # sigma unreachable
        else:
            if distance_from_shell < threshold_distance:
                # near a;
                if __debug__:
                    log.debug('GF: only a')
                return GreensFunction3DAbs(self.D_r, r0, self.a_r)
                
            else:
                # distant from both a and sigma; 
                if __debug__:
                    log.debug('GF: free')
                return GreensFunction3D(self.D_r, r0)

    def create_com_vector(self, r):
        return random_vector(r)

    def create_interparticle_vector(self, gf, r, dt, r0, old_iv): 

        if __debug__:
            log.debug("create_interparticle_vector: r=%g, dt=%g", r, dt)
        theta = draw_theta_wrapper(gf, r, dt)

        if theta == 0.0:
            # no rotation is necessary -> new_iv = new_iv
            new_iv = old_iv
        elif theta % numpy.pi != 0.0:
            # alternative calculation
            # rotate the old_iv to the new theta
            rotation_axis = crossproduct_against_z_axis(old_iv)
            rotation_axis = normalize(rotation_axis)
            new_iv = rotate_vector(old_iv, rotation_axis, theta)

            # rotate the new_iv around the old_iv with angle phi
            phi = myrandom.uniform() * 2 * Pi
            rotation_axis = normalize(old_iv)
            new_iv = rotate_vector(new_iv, rotation_axis, phi)
        else:
            # theta == numpi.pi -> just mirror the old_iv 
            new_iv = -old_iv

        # adjust length of the vector
        # note that r0 = length (old_iv)
        new_iv = (r/r0) * new_iv

        return new_iv

    def __str__(self):
        return 'Spherical' + Pair.__str__(self)


#class PlanarSurfacePair(SimplePair):
class PlanarSurfacePair(SimplePair, Others, hasCylindricalShell):
    """2 Particles inside a (cylindrical) shell on a PlanarSurface. 
    (Hockey pucks).

    """
    def __init__(self, domain_id, shell_id, testShell, rrs):
#    def __init__(self, domain_id, shell_center, single1, single2, shell_id,
#                 r0, shell_size, rrs, structure):
#        single1 = testShell.single1
#        single2 = testShell.single2
#        structure = testShell.structure
#        shell_center = testShell.center
#        shell_size = testShell.radius  # TODO this should be removed
#        r0 = testShell.r0        # TODO this should be removed

#        SimplePair.__init__(self, domain_id, shell_center, single1, single2, shell_id,
#                      r0, shell_size, rrs, structure)
        hasCylindricalShell.__init__(self, testShell, domain_id)
        SimplePair.__init__(self, domain_id, shell_id, rrs)

    def com_greens_function(self):
        return GreensFunction2DAbsSym(self.D_R, self.a_R)

    def iv_greens_function(self, r0):
        return GreensFunction2DRadAbs(self.D_r, self.interparticle_ktot, r0,
                                      self.sigma, self.a_r)

#    def create_new_shell(self, position, radius, domain_id):
#        # The half_length (thickness/2) of a hockey puck is not more 
#        # than it has to be (namely the radius of the particle), so if 
#        # the particle undergoes an unbinding reaction we still have to 
#        # clear the target volume and the move may be rejected (NoSpace 
#        # error).
#        orientation = crossproduct(self.structure.shape.unit_x,
#                                   self.structure.shape.unit_y)
#        half_length = max(self.pid_particle_pair1[1].radius,
#                          self.pid_particle_pair2[1].radius)
#        return CylindricalShell(domain_id, Cylinder(position, radius, 
#                                                    orientation, half_length))

    def get_shell_direction(self):
        return self.shell.shape.unit_z

    def choose_pair_greens_function(self, r0, t):
        # selects between the full solution or an approximation where one of
        # the boundaries is ignored
        distance_from_sigma = r0 - self.sigma
        distance_from_shell = self.a_r - r0

        threshold_distance = Pair.CUTOFF_FACTOR * \
            math.sqrt(4.0 * self.D_r * t)

        if distance_from_sigma < threshold_distance:
        
            if distance_from_shell < threshold_distance:
                # near both a and sigma;
                # use GreensFunction2DRadAbs
                if __debug__:
                    log.debug('GF2D: normal')
                return self.iv_greens_function(r0)
            else:
                # near sigma; use GreensFunction2DRadInf
                if __debug__:
                    log.debug('GF2D: only sigma')
                return self.iv_greens_function(r0)
                # Todo
                # return GreensFunction2DRadInf(self.D_r, self.interparticle_ktot, r0, self.sigma)
        else:
            if distance_from_shell < threshold_distance:
                # near a;
                if __debug__:
                    log.debug('GF2D: only a')
                return self.iv_greens_function(r0)
                # Todo
                # return GreensFunction2DAbs(self.D_r, r0, self.a_r)
                
            else:
                # distant from both a and sigma; 
                if __debug__:
                    log.debug('GF2D: free')
                return self.iv_greens_function(r0)
                # Todo
                # return GreensFunction2D(self.D_r, r0)

    def create_com_vector(self, r):
        x, y = random_vector2D(r)
        # project the com onto the surface unit vectors to make sure that the coordinates are in the surface
        return x * self.structure.shape.unit_x + y * self.structure.shape.unit_y

    def create_interparticle_vector(self, gf, r, dt, r0, old_iv): 

        if __debug__:
            log.debug("create_interparticle_vector: r=%g, dt=%g", r, dt)
        theta = draw_theta_wrapper(gf, r, dt)

        # note that r0 = length (old_iv)
        new_iv = (r/r0) * rotate_vector(old_iv, self.structure.shape.unit_z, theta)
        # note that unit_z can point two ways rotating the vector clockwise or counterclockwise
        # Since theta is symmetric this doesn't matter.

        # project the new_iv down on the unit vectors of the surface to prevent the particle from
        # leaving the surface due to numerical problem
        unit_x = self.structure.shape.unit_x
        unit_y = self.structure.shape.unit_y

        new_iv_x = unit_x * numpy.dot(new_iv, unit_x)
        new_iv_y = unit_y * numpy.dot(new_iv, unit_y)

        return new_iv_x + new_iv_y

    def __str__(self):
        return 'PlanarSurface' + Pair.__str__(self)


#class CylindricalSurfacePair(SimplePair):
class CylindricalSurfacePair(SimplePair, Others, hasCylindricalShell):
    """2 Particles inside a (cylindrical) shell on a CylindricalSurface.  
    (Rods).

    """
    def __init__(self, domain_id, shell_id, testShell, rrs):
#    def __init__(self, domain_id, shell_center, single1, single2, shell_id,
#                 r0, shell_size, rrs, structure):
#        SimplePair.__init__(self, domain_id, shell_center, single1, single2, shell_id,
#                      r0, shell_size, rrs, structure)
        hasCylindricalShell.__init__(self, testShell, domain_id)
        SimplePair.__init__(self, domain_id, shell_id, rrs)

    def get_v_tot(self):
        return self.pid_particle_pair2[1].v - \
               self.pid_particle_pair1[1].v
    v_tot = property(get_v_tot)
    
    v_r = property(get_v_tot)

    def get_v_R(self):
        return (self.pid_particle_pair1[1].v * 
                self.pid_particle_pair2[1].D +
                self.pid_particle_pair2[1].v *
                self.pid_particle_pair1[1].D) / self.D_tot
    v_R = property(get_v_R)

    def com_greens_function(self):
        # The domain is created around r0, so r0 corresponds to r=0 within the domain
        return GreensFunction1DAbsAbs(self.D_R, self.v_R, 0.0, -self.a_R, self.a_R)

    def iv_greens_function(self, r0):
        return GreensFunction1DRadAbs(self.D_r, self.v_r, self.interparticle_ktot, r0, self.sigma, self.a_r)

#    def create_new_shell(self, position, half_length, domain_id):
#        # The radius of a rod is not more than it has to be (namely the 
#        # radius of the biggest particle), so if the particle undergoes 
#        # an unbinding reaction we still have to clear the target volume 
#        # and the move may be rejected (NoSpace error).
#        radius = max(self.pid_particle_pair1[1].radius,
#                     self.pid_particle_pair2[1].radius)
#        orientation = self.structure.shape.unit_z
#        return CylindricalShell(domain_id, Cylinder(position, radius, 
#                                                    orientation, half_length))

    def get_shell_direction(self):
        return self.shell.shape.unit_z

    def choose_pair_greens_function(self, r0, t):
        # Todo
        return self.iv_greens_function(r0)

    def create_com_vector(self, r):
        # project the com onto the surface unit vector to make sure that the coordinates are in the surface
        return r * self.structure.shape.unit_z

    def create_interparticle_vector(self, gf, r, dt, r0, old_iv): 
        if __debug__:
            log.debug("create_interparticle_vector: r=%g, dt=%g", r, dt)

        # note that r0 = length (old_iv)
        new_iv = (r/r0) * old_iv

        # project the new_iv down on the unit vectors of the surface to prevent the particle from
        # leaving the surface due to numerical problem
        unit_z = self.structure.shape.unit_z
        new_iv_z = unit_z * numpy.dot(new_iv, unit_z)

        return new_iv_z

#    def get_shell_size(self):
#        # Heads up.
#        return self.shell.shape.half_length

    def __str__(self):
        return 'CylindricalSurface' + Pair.__str__(self)


class MixedPair2D3D(Pair, Others, hasCylindricalShell):

#    @classmethod
#    def get_min_shell_dimensions(cls, single1, single2, geometrycontainer):
#
#        assert isinstance (single1.structure, PlanarSurface)
#        assert isinstance (single2.structure, CuboidalRegion)
#
#        pos1 = single1.pid_particle_pair[1].position
#        pos2 = single2.pid_particle_pair[1].position
#        radius1 = single1.pid_particle_pair[1].radius
#        radius2 = single2.pid_particle_pair[1].radius
#        D1 = single1.pid_particle_pair[1].D
#        D2 = single2.pid_particle_pair[1].D
#        D12 = D1 + D2
#        z_scaling = cls.calc_z_scaling_factor (D1, D2)
#
#        # calculate the interparticle vector and project it onto the coordinate system of the structure
#        # of particle1
#        pos2_t = geometrycontainer.world.cyclic_transpose(pos2, pos1)
#        iv = pos2_t - pos1
#        iv_x =     numpy.dot(iv, single1.structure.shape.unit_x)
#        iv_y =     numpy.dot(iv, single1.structure.shape.unit_y)
#        iv_z = abs(numpy.dot(iv, single1.structure.shape.unit_z))
#
#        r0 = math.sqrt( iv_x**2 + iv_y**2 + (iv_z*z_scaling)**2)
#
#        # calculate the minimal height z_right1 of the shell including burst radius
#        # with the accompanying radius r1
#        z_right1 = iv_z + radius2 * SINGLE_SHELL_FACTOR
#        r1       = cls.r_right(single1, single2, r0, z_right1)
#
#        len_iv = r0
#        iv_shell_radius1 = len_iv * D1 / D12 
#        iv_shell_radius2 = len_iv * D2 / D12 
#
#        # fix the minimum shell size for the CoM domain
#        com_shell_radius = max(radius1, radius2)
#
#        # calculate the minimal dimensions of the protective domain including space for the
#        # burst volumes of the particles
#        r2 = max(iv_shell_radius1 + com_shell_radius + radius1 * SINGLE_SHELL_FACTOR,
#                 iv_shell_radius2 + com_shell_radius + radius2 * SINGLE_SHELL_FACTOR)
#        z_right2 = cls.z_right(single1, single2, r0, r2)
#
#        # of both alternatives pick the largest one
#        z_right, r = max((z_right1, r1),(z_right2, r2))
#
#        # z_left is always the same
#        z_left = radius1
#
#        return r, z_left, z_right
#
#    @classmethod
#    def get_max_shell_dimensions(cls, single1, single2, geometrycontainer, domains):
#        # Properties of the MixedPair system
#
#        assert isinstance (single1.structure, PlanarSurface)
#        assert isinstance (single2.structure, CuboidalRegion)
#
#        com, iv = cls.do_transform(single1, single2, geometrycontainer.world)
#        reference_point = com
#        orientation_vector = normalize(single1.structure.shape.unit_z * numpy.dot (iv, single1.structure.shape.unit_z))
#        r0 = length(iv)
#
#        # Make sure the maximal cylinder fits in the maximal sphere. Matrix space
#        # doesn't allow to check for shells outside the maximal sphere.
#        max_cylinder_radius      = geometrycontainer.get_max_shell_size()/math.sqrt(2)
#        max_cylinder_half_length = max_cylinder_radius
#
#        # set initial values of shell dimensions
#        z_left = single1.pid_particle_pair[1].radius
#        z_right = max_cylinder_half_length * 2 - z_left
#        r = max_cylinder_radius
#        
#        # get all the domains in the neighborhood
#        search_point = reference_point + ((z_right - z_left)/2.0) * orientation_vector
#        search_point = geometrycontainer.world.apply_boundary(search_point)
#        all_neighbor_ids = \
#            geometrycontainer.get_neighbors_within_radius_no_sort(search_point,
#                                                     geometrycontainer.get_max_shell_size(),
#                                                     ignore=[single1.domain_id, single2.domain_id])
#        all_neighbors = [domains[domain_id] for domain_id in all_neighbor_ids]
#
#        world = geometrycontainer.world
#
#        # adjust the maximal shell dimensions for every shell in the neighborhood
#        for domain in all_neighbors:
#            if isinstance(domain, Multi):
#                for _, shell in domain.shell_list:
#                    shell_position = shell.shape.position
#                    shell_size = shell.shape.radius
#                    r, z_left, z_right = cls.laurens_algorithm(single1, single2, r0, shell_position, shell_size, shell,
#                                                                     reference_point, orientation_vector, r, z_left,
#                                                                     z_right, world)
#            else:
#                shell_position = domain.shell.shape.position
#                shell_size = domain.get_shell_size()
#
#                if domain.dt == 0.0 and domain.getD() > 0:
#                    # This is one of the bursted singles.
#                    # Or a particle that just escaped it's multi.
#                    shell_size *= SINGLE_SHELL_FACTOR
#
#                r, z_left, z_right = cls.laurens_algorithm(single1, single2, r0, shell_position, shell_size, domain.shell,
#                                                                 reference_point, orientation_vector, r, z_left,
#                                                                 z_right, world)
#
#
#
#        return r, z_left, z_right

    ################
    @classmethod
    def calc_z_scaling_factor(cls, D2d, D3d):
        # calculates the scaling factor to make the anisotropic diffusion problem into a isotropic one
        return math.sqrt( (D2d + D3d) / D3d)

    def __init__(self, domain_id, shell_id, testShell, rrs):
#    def __init__(self, domain_id, single1, single2, shell_id, shell_center, shell_radius,
#                 shell_half_length, shell_orientation_vector, r0, rrs):

#        Pair.__init__(self, domain_id, single1, single2, shell_id, r0, rrs)

        hasCylindricalShell.__init__(self, testShell, domain_id)
        Pair.__init__(self, domain_id, shell_id, rrs)

        assert self.testShell
        assert isinstance(self.testShell, testMixedPair2D3D)
        self.surface = self.testShell.surface       # the surface on which particle1 lives
        self.structure = self.testShell.structure   # structure on which both particles live

        self.z_scaling_factor = self.testShell.get_scaling_factor()

        # set the radii of the inner domains as a function of the outer protective domain
#        self.r0 = r0
        self.a_R, self.a_r = self.determine_radii()

#        self.shell = self.create_new_shell (shell_center, shell_radius, shell_half_length, shell_orientation_vector,
#                                            domain_id)

    def determine_radii(self):
        # determines the dimensions of the domains used for the Green's functions from the dimensions
        # of the cylindrical shell.
        # Note that the function assumes that the cylinder is dimensioned properly
        radius1 = self.pid_particle_pair1[1].radius
        radius2 = self.pid_particle_pair2[1].radius
        D_1 = self.pid_particle_pair1[1].D
        D_2 = self.pid_particle_pair2[1].D

        half_length = self.shell.shape.half_length
        radius = self.shell.shape.radius

        # use class methods to check dimensions of the cylinder
        # note that we assume that only the z_right is scalable
        r_check   = self.testShell.r_right (half_length*2 - radius1)
        hl2_check = self.testShell.z_right (radius) + self.testShell.z_left  (radius)
        if r_check > radius:
            # radius should have been larger, adjust z_right instead
            if __debug__:
                log.debug('MixedPair2D3D: half_length was too high for radius. '
                          'radius = %.3g, r_check = %.3g, 2half_length = %.3g, 2hl_check = %.3g' %
                          (radius, r_check, half_length*2, hl2_check))

            half_length = hl2_check/2

        elif hl2_check > (half_length*2):
            # half_length should have been larger, adjust radius instead
            if __debug__:
                log.debug('MixedPair2D3D: radius was too high for half_length. '
                          'radius = %.3g, r_check = %.3g, 2half_length = %.3g, 2hl_check = %.3g' %
                          (radius, r_check, half_length*2, hl2_check))
            radius = r_check
        else:
            # dimensions were ok
            pass

        z_left  = radius1
        z_right = half_length * 2 - z_left
        a_r = (z_right - radius2) * self.z_scaling_factor

        # calculate the maximal displacement of a particle given an a_r. Also include its radius
        space_for_iv = max( (D_1/self.D_tot) * a_r + radius1,
                            (D_2/self.D_tot) * a_r + radius2)
        a_R = radius - space_for_iv


        if __debug__:
            tr = ((a_r - self.r0)**2) / (6 * self.D_r)  # the expected escape time of the iv
            if self.D_R == 0:
                tR = numpy.inf 
            else:
                tR = (a_R**2) / (4 * self.D_R)          # the expected escape time of the CoM
            log.debug('a_r= %s, tr= %s, a_R= %s, tR= %s' % \
                      (FORMAT_DOUBLE % a_r, FORMAT_DOUBLE % tr,
                       FORMAT_DOUBLE % a_R, FORMAT_DOUBLE % tR))

        assert (self.sigma < a_r) and (a_r < half_length*2)
        assert (0 < a_R) and (a_R < radius)

        return a_R, a_r

#    def get_com(self):
#        # Note that it is implied that the center of the shell is chosen such that its projection on the
#        # membrane is the CoM of the particles.
#        dist_to_surface = self.shell.shape.half_length - self.pid_particle_pair1[1].radius
#        return self.shell.shape.position - dist_to_surface * self.shell.shape.unit_z
#    com = property(get_com)
#
#    def get_scaling_factor(self):
#        D_2 = self.pid_particle_pair2[1].D      # particle 2 is in 3D and is the only contributor to diffusion
#                                                # normal to the plane
#        return math.sqrt(self.D_r/D_2)
#    z_scaling_factor = property(get_scaling_factor)
#
#    @classmethod
#    def do_transform(cls, single1, single2, world):
#
#        pos1 = single1.pid_particle_pair[1].position
#        pos2 = single2.pid_particle_pair[1].position
#        D_1 = single1.pid_particle_pair[1].D
#        D_2 = single2.pid_particle_pair[1].D
#
#        surface = single1.structure
#        # we assume that the particle of single1 lives on the membrane
#        assert isinstance(surface, PlanarSurface)
#
#
#        # the CoM is calculated in a similar way to a normal 3D pair
#        com = (D_2 * pos1 + D_1 * pos2) / (D_1 + D_2)
#        # and then projected onto the plane to make sure the CoM is in the surface
#        com = world.cyclic_transpose(com, surface.shape.position)
#        com, _ = surface.projected_point (com)
#        com = world.apply_boundary(com)
#
#        # calculate the interparticle vector
#        pos2t = world.cyclic_transpose(pos2, pos1)
#        iv = pos2t - pos1
#
#        # calculate the scaling factor due to anisotropic diffusion of iv
#        z_scaling_factor = math.sqrt((D_1 + D_2)/D_2)
#
#        # move and recale the iv in the axis normal to the plane
#        iv_z = surface.shape.unit_z * numpy.dot (iv, surface.shape.unit_z)
#        iv_z_length = length(iv_z)
#
#        new_iv_z_length = (iv_z_length) * z_scaling_factor
#        new_iv_z = (new_iv_z_length / iv_z_length) * iv_z
#
#        iv = iv - iv_z + new_iv_z
#
#        return com, iv

    @ classmethod
    def do_back_transform(cls, com, iv, D1, D2, radius1, radius2, surface, unit_z):
    # here we assume that the com and iv are really in the structure and no adjustments have to be
    # made

        D_tot = D1 + D2
        weight1 = D1 / D_tot
        weight2 = D2 / D_tot

#        min_iv_z_length = radius2
        # get a sence of scale of the separation between the particle and the membrane
        min_iv_z_length = (radius1 + radius2) * 0.5 * (MINIMAL_SEPARATION_FACTOR - 1.0)

        # get the coordinates of the iv relative to the system of the surface (or actually the shell)
        iv_x = surface.shape.unit_x * numpy.dot(iv, surface.shape.unit_x)
        iv_y = surface.shape.unit_y * numpy.dot(iv, surface.shape.unit_y)

        # reflect the coordinates in the unit_z direction back to the side of the membrane
        # where the domain is. Note that it's implied that the origin of the coordinate system lies in the
        # plane of the membrane
        iv_z_length = abs(numpy.dot(iv, unit_z))   # FIXME maybe first project the shell unit_z onto the 
                                                   # surface unit_z to prevent numerical problems?
        # do the reverse scaling
        iv_z_length = iv_z_length / math.sqrt(D_tot/D2)

        # if the particle is overlapping with the membrane, make sure it doesn't
        # IS THIS STILL NECESSARY?
        if iv_z_length < min_iv_z_length:
            iv_z_length = min_iv_z_length #* MINIMAL_SEPARATION_FACTOR

        iv_z = unit_z * iv_z_length

        pos1 = com - weight1 * (iv_x + iv_y)
        pos2 = com + weight2 * (iv_x + iv_y) + iv_z 

        return pos1, pos2

#    def get_sigma(self):
#        # rescale sigma to correct for the rescaling of the coordinate system
#        # This is the sigma that is used for the evaluation of the Green's function and is in this case slightly
#        # different than the sums of the radii of the particleso
#        xi = self.z_scaling_factor
#        xi_inv = 1.0/xi
#        alpha = math.acos(xi_inv)
#        sigma = self.pid_particle_pair1[1].radius + self.pid_particle_pair2[1].radius
#        rho = abs(sigma * math.sqrt(0.5 + (alpha * xi/(2.0*math.sin(alpha)))))
#        return rho
#    sigma = property (get_sigma)

    def choose_pair_greens_function(self, r0, t):
        return self.iv_greens_function(r0)

    def com_greens_function(self):
        # diffusion of the CoM of a mixed pair is only two dimensional
        return GreensFunction2DAbsSym(self.D_R, self.a_R)

    def iv_greens_function(self, r0):
        # diffusion of the interparticle vector is three dimensional
        # TODO Fix ugly hack to prevent particle overlap problem
        return GreensFunction3DRadAbs(self.D_r, self.interparticle_ktot, max(r0, self.sigma),
                                      self.sigma, self.a_r)

#    def create_new_shell(self, position, radius, half_length, orientation_vector, domain_id):
#        # It is assumed that the parameter are correct -> the 'bottom' of the cylinder
#        # sticks through the membrane enough to accomodate the particle in the membrane
#        # assert that orientation_vector is on same line as surface.shape.unit_z
#        return CylindricalShell(domain_id, Cylinder(position, radius, orientation_vector,
#                                                    half_length))
#
#    def get_shell_size(self):
#        return self.shell.shape.radius

    def get_shell_direction(self):
        return self.shell.shape.unit_z

    def create_com_vector(self, r):
        x, y = random_vector2D(r)
        # project the com onto the surface unit vectors to make sure that the coordinates are in the surface
        return x * self.surface.shape.unit_x + y * self.surface.shape.unit_y

    def create_interparticle_vector(self, gf, r, dt, r0, old_iv):
        # FIXME This is actually the same method as for SphericalPair

        if __debug__:
            log.debug("create_interparticle_vector: r=%g, dt=%g", r, dt)

        theta = draw_theta_wrapper(gf, r, dt)   # draw_theta return a value between -Pi and Pi

        if theta == 0.0:
            new_iv = old_iv
        elif theta % numpy.pi != 0.0:
            # rotate the old_iv to the new theta
            rotation_axis = crossproduct_against_z_axis(old_iv)
            rotation_axis = normalize(rotation_axis)
            new_iv = rotate_vector(old_iv, rotation_axis, theta)

            # rotate the new_iv around the old_iv with angle phi
            phi = myrandom.uniform() * 2 * Pi
            rotation_axis = normalize(old_iv)
            new_iv = rotate_vector(new_iv, rotation_axis, phi)
        else:
            # theta == Pi -> just mirror the old_iv
            new_iv = -old_iv

        # adjust length of the vector to new r
        # note that r0 = length (old_iv)
        new_iv = (r/r0) * new_iv

        return new_iv

    def __str__(self):
        return 'Mixed' + Pair.__str__(self)

class MixedPair1D3D(Pair):

    @classmethod
    def do_back_transform(cls, com, iv, D1, D2, radius1, radius2, surface, unit_z):
    # here we assume that the com and iv are really in the structure and no adjustments have to be
    # made
    # unit_z is a putatively different unit_z than that available in the surface object

        D_tot = D1 + D2
        weight1 = D1 / D_tot
        weight2 = D2 / D_tot

        min_iv_r_length = radius2 + surface.shape.radius

        # get the coordinates of the iv relative to the system of the surface (or actually the shell)
        iv_z = surface.shape.unit_z * numpy.dot(iv, surface.shape.unit_z)
        iv_r = iv - iv_z

        # reflect the coordinates in the unit_z direction back to the side of the membrane
        # where the domain is. Note that it's implied that the origin of the coordinate system lies in the
        # plane of the membrane
        iv_r_length = length(iv_r)
        # do the reverse scaling
        iv_r_length_new = iv_r_length / cls.calc_r_scaling_factor(D1, D2)

        # if the particle is overlapping with the cylinder, make sure it doesn't
        if iv_r_length_new < min_iv_r_length:
            iv_r_length_new = min_iv_r_length * MINIMAL_SEPARATION_FACTOR

        iv_r =  (iv_r_length_new/iv_r_length) * iv_r

        pos1 = com - weight1 * iv_z
        pos2 = com + weight2 * iv_z + iv_r

        return pos1, pos2

    @classmethod
    def calc_r_scaling_factor(cls, D1, D2):
        return math.sqrt((D1 + D2)/D2)

    def get_shell_direction(self):
        return self.shell.shape.unit_z

