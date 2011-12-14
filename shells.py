#!/usr/bin/python


class Others(object):

    # stuff for hasShell subclasses
    def shell_list_for_single(self):
        return self.shell_list
    def shell_list_for_other(self):
        return self.shell_list

    # stuff for testShell subclasses
    def get_neighbor_shell_list(self, neighbor):
        # note that this returns the list of shells, NOT dr, dz_right, dz_left
        return neighbor.shell_list_for_other()

class NonInteractionSingles(object):

    # stuff for hasShell subclasses
    def shell_list_for_single(self):
        pass    # needs to be overloaded
    def shell_list_for_other(self):
        pass    # needs to be overloaded

    # stuff for testShell subclasses
    def get_neighbor_shell_list(self, neighbor):
        # note that this returns the list of shells, NOT dr, dz_right, dz_left
        return neighbor.shell_list_for_single()

##################################
class hasShell(object):

    def shell_list(self):
            return self.shell_list

#######################
class hasSphericalShell(hasShell):

    self.shell

    def __init__(self, sphericaltestShell):
        self.center = sphericaltestShell.center
        self.radius = sphericaltestShell.radius

    # methods to size up a shell to this shell
    def get_dr_dzright_dzleft_to_shell(self, shell, domain_to_scale, r, z_right, z_left):
        # This will scale the 'domain_to_scale' (a cylinder) using the 'shell_appearance' as the limiting shell
        # CylindricaltestShell ('domain_to_scale') -> SphericalShell ('self' with 'shell_appearance')

        assert type(shell, Sphere)        # because the shell actually is part of the current object

        # Do Laurens' algorithm (part1)
        shell_position = shell.shape.position
        shell_size = shell.shape.radius

        # get the reference point and orientation of the domain to scale
        reference_point = domain_to_scale.get_referencepoint()
        orientation_vector = domain_to_scale.get_orientation_vector()

        # determine on what side the midpoint of the shell is relative to the reference_point
        shell_position_t = world.cyclic_transpose (shell_position, reference_point)
        ref_to_shell = shell_position_t - reference_point
        ref_to_shell_z_len = numpy.dot(ref_to_shell, orientation_vector)
        # express the ref_to_shell vector in the coordinate system
        ref_to_shell_z = ref_to_shell_z_len * orientation_vector
        ref_to_shell_r = ref_to_shell - ref_to_shell_z

        # if the shell is on the side of orientation_vector -> use z_right
        # also use this one if the shell is in plane with the reference_point
        if ref_to_shell_z_len >= 0:
            direction = 1           # improve direction specification such that following calculations still work
                                    # Otherwise problems when direction = 0
            get_scale_angle  = domain_to_scale.get_right_scalingangle
            get_scale_center = domain_to_scale.get_right_scalingcenter
            r1_function =  cls.r_right
            z1_function = cls.z_right
            z2_function = cls.z_left
            z1 = z_right
            z2 = z_left
        # else -> use z_left
        else:
            direction = -1
            get_scalingangle  = domain_to_scale.get_left_scalingangle
            get_scalingcenter = domain_to_scale.get_left_scalingcenter
            r1_function =  cls.r_left
            z1_function = cls.z_left
            z2_function = cls.z_right
            z1 = z_left
            z2 = z_right

        # correct the orientation vectors
        orientation_z = orientation_vector * direction
        if length(ref_to_shell_r) == 0:
            # produce a random vector perpendicular to the 'orientation_vector_z'
            unit_vector3D = random_unit_vector()    # TODO find cheaper way
            orientation_r = normalize(unit_vector3D - numpy.dot(unit_vector3D, orientation_z))
        else:
            orientation_r = normalize(ref_to_shell_r)


        # calculate the center from which linear scaling will take place
        phi = get_scalingangle()
        scale_center = get_scalingcenter(orientation_z, orientation_r) + reference_point
        scale_center = world.apply_boundary(scale_center)


        # calculate the vector from the scale center to the center of the shell
        shell_position_t = world.cyclic_transpose (shell_position, scale_center)
        shell_scale_center = shell_position_t - scale_center
        shell_scalecenter_z = numpy.dot(shell_scale_center, orientation_vector_z)
        shell_scalecenter_r = numpy.dot(shell_scale_center, orientation_vector_r)


        # calculate the angle theta of the vector from the scale center to the shell with the vector
        # to the scale center (which is +- the orientation_vector)
        theta = vector_angle (shell_scale_center, orientation_z)

        psi = theta - phi

        ### check what situation arrises
        # The shell can hit the cylinder on the flat side (situation 1),
        #                                on the edge (situation 2),
        #                                or on the round side (situation 3).
        # I think this also works for phi == 0 and phi == Pi/2
        if psi <= -phi:
        # The midpoint of the shell lies on the axis of the cylinder
            situation = 1

        elif -phi < psi and psi < 0:
        # The (spherical) shell can touch the cylinder on its flat side or its edge
            if phi == Pi/2.0:
                r_tan_phi = 0
            else:
                r_tan_phi = abs(shell_scalecenter_r)/math.tan(phi)   # phi == 0 should not get here
                                                                     # TODO the abs may be wrong/unnecessary

            a_thres = shell_scalecenter_z - r_tan_phi
            if shell_size < a_thres:
                situation = 1
            else:
                situation = 2

        elif 0 <= psi and psi < (Pi/2.0 - phi):
        # The (spherical) shell can touch the cylinder on its edge or its radial side
            tan_phi = math.tan(phi)                             # phi == Pi/2 should not get here
            a_thres = abs(shell_scalecenter_r) - shell_scalecenter_z * tan_phi  # TODO same here
            if shell_size > a_thres:
                situation = 2
            else:
                situation = 3

        elif (Pi/2.0 - phi) <= psi:
        # The shell is always only on the radial side of the cylinder
            situation = 3

        else:
        # Don't know how we would get here, but it shouldn't happen
            raise RuntimeError('Error: psi was not in valid range. psi = %s, phi = %s, theta = %s' %
                               (FORMAT_DOUBLE % psi, FORMAT_DOUBLE % phi, FORMAT_DOUBLE % theta))


        ### Get the right values for z and r for the given situation
        if situation == 1:      # shell hits cylinder on the flat side
            z1_new = min(z1, (scalecenter_h0 + shell_scalecenter_z - shell_size)/SAFETY)
            r_new  = min(r,  r1_function(single1, single2, r0, z1_new))
            z2_new = min(z2, z2_function(single1, single2, r0, r_new))

        elif situation == 2:    # shell hits sphere on the edge
            a_sq = shell_size*shell_size
            shell_scalecenter_len = length(shell_scale_center)
            ss_sq = shell_scalecenter_len*shell_scalecenter_len
            sin_phi = math.sin(phi)
            sin_psi = math.sin(psi)
            cos_psi = math.cos(psi)
            scalecenter_shell_dist = (shell_scalecenter_len * cos_psi - math.sqrt(a_sq - ss_sq*sin_psi*sin_psi) )
            # FIXME UGLY FIX BELOW
            scalecenter_shell_dist /= 1.1
            scalecenter_shell_dist /= SAFETY

            if phi <= Pi/4:
                z1_new = min(z1, scalecenter_h0 + cos_phi * scalecenter_shell_dist)
                r_new  = min(r, r1_function(single1, single2, r0, z1_new))
                z2_new = min(z2, z2_function(single1, single2, r0, r_new))
            else:
                r_new = min(r, scalecenter_r0 + sin_phi * scalecenter_shell_dist)
                z1_new = min(z1, z1_function(single1, single2, r0, r_new))
                z2_new = min(z2, z2_function(single1, single2, r0, r_new))

        elif situation == 3:    # shell hits cylinder on the round side
            r_new = min(r, (scalecenter_r0 + abs(shell_scalecenter_r) - shell_size)/SAFETY)
            z1_new = min(z1, z1_function(single1, single2, r0, r_new))
            z2_new = min(z2, z2_function(single1, single2, r0, r_new))
        else:
            raise RuntimeError('Bad situation for MixedPair shell making')

        # switch the z values in case it's necessary. r doesn't have to be switched.
        r = r_new
        if direction >= 0.0:
            z_right = z1_new
            z_left  = z2_new
        else:
            z_right = z2_new
            z_left  = z1_new

        return (r, z_right, z_left)


    def get_radius_to_shell(self, shell_appearance, domain_to_scale):
        # SphericaltestShell ('domain_to_scale') -> SphericalShell ('self' with 'shell_appearance')

        assert type(shell_appearance, Sphere)        # because the shell actually originated here
        # Do simple distance calculation to sphere

        return r

#####
class SphericalSingle(hasSphericalShell, NonInteractionSingles):

    # these return potentially corrected dimensions
    def shell_list_for_single(self):
        min_radius = self.get_min_radius()
        if  self.shell.shape.radius < min_radius:
            fake_shell = make_new_shell (min_radius)
            return [fake_shell, ]
        else:
            return self.shell_list

    def shell_list_for_other(self):
        min_radius = particle_radius * SINGLE_SHELL_FACTOR
        if self.shell.shape.radius < min_radius:
            fake_shell = make_new_shell (min_radius)
            return [fake_shell, ]
        else:
            return self.shell_list

#####
class SphericalPair(hasSphericalShell, Others):

#####
class Multi(hasSphericalShell, Others):


#########################
class hasCylindricalShell(hasShell):

    def __init__(self, cylindricaltestShell):
        self.center = uitrekenen zie interactions
        self.radius = cylindricaltestShell.dr
        self.half_length = uitrekenen zie interactions

    def get_center:

    def get_dr_dzright_dzleft_to_shell(self, shell_appearance, domain_to_scale):
        # This will scale the 'domain_to_scale' (a cylinder) using the 'shell_appearance' as the limiting shell
        # CylindricaltestShell ('domain_to_scale') -> CylindricalShell ('self' with 'shell_appearance')

        assert type(shell_appearance, Cylinder)        # because the shell actually originated here
        # Do Laurens' algorithm (part2)

        return (r, z_right, z_left)

    def get_radius_to_shell(self, shell_appearance, domain_to_scale):
        # SphericaltestShell ('domain_to_scale') -> CylindricalShell ('self' with 'shell_appearance')

        assert type(shell_appearance, Cylinder)        # because the shell actually originated here
        # Do simple distance calculation to cylinder

        return r

#####
class PlanarSurfaceSingle(hasCylindricalShell, NonInteractionSingles):

    # these return corrected dimensions, since we reserve more space for the NonInteractionSingle
    def shell_list_for_single(self):
        min_radius = particle_radius * MULTI_SHELL_FACTOR
        if self.shell.shape.radius < min_radius:
            fake_shell = create_shell(min_radius, normal_half_length)
            return [fake_shell, ]
        else:
            return self.shell_list

    def shell_list_for_other(self):
        min_radius = particle_radius * SINGLE_SHELL_FACTOR
        if self.shell.shape.radius < min_radius:
            fake_shell = create_shell(min_radius, normal_half_length)
            return [fake_shell, ]
        else:
            return self.shell_list

#####
class CylindricalSurfaceSingle(hasCylindricalShell, NonInteractionSingles):

    # these return corrected dimensions, since we reserve more space for the NonInteractionSingle
    def shell_list_for_single(self):
        min_half_length = particle_radius * MULTI_SHELL_FACTOR
        if self.shell.shape.half_length < min_half_length:
            fake_shell = create_shell(radius, min_half_length)
            return [fake_shell, ]
        else:
            return self.shell_list

    def shell_list_for_other(self):
        min_half_length = particle_radius * SINGLE_SHELL_FACTOR
        if self.shell.shape.half_length < min_half_length:
            fake_shell = create_shell(radius, min_half_length)
            return [fake_shell, ]
        else:
            return self.shell_list

#####
class PlanarSurfacePair(hasCylindricalShell, Others):
class PlanarSurfaceInteraction(hasCylindricalShell, Others):
class CylindricalSurfacePair(hasCylindricalShell, Others):
class CylindricalSurfaceInteraction(hasCylindricalShell, Others):
class MixedPair3D2D(hasCylindricalShell, Others):
#class MixedPair3D1D(hasCylindricalShell):


#######################################################
#######################################################
class testShell(object):

    def get_neighbors(self):
        searchpoint = self.get_searchpoint()
        radius = self.get_searchradius()
        return self.geometrycontainer.get_neigbors(searchpoint, radius)

########################
class SphericaltestShell(testShell):

    self.radius
    self.center

    def determine_possible_shell(self):
        neighbors = self.get_neighbors()
        min_radius = self.get_min_radius()
        max_radius = self.get_max_radius()

        radius = max_radius
        for neighbor in neighbors:

            shell_list = self.get_neighbor_shell_list(neighbor)
            for _, shell_appearance in shell_list:
                new_radius = neighbor.get_radius_to_shell(shell_appearance, self)
                # if the newly calculated dimensions are smaller than the current one, use them
                radius = min(radius, new_radius)

                if radius < min_radius:
                    radius = None
                    break

        if radius:
            return self.create_new_shell(radius)
        else
            return None

    def get_searchpoint(self):
        return self.center
    def get_searchradius(self):
        return self.geometrycontainer.max()

    def get_max_radius(self):
        return shell_container.get_max_shell_size()

#####
class SphericalSingletestShell(testSphericalShell, NonInteractionSingles):

    self.center = particle_position

    def get_min_radius(self):
        return particle_radius * MULTI_SHELL_FACTOR     # the minimum radius of a NonInteractionSingle

#####
class SphericalPairtestShell(testSphericalShell, Others):

    self.center = calculate CoM
    def get_min_radius:
        return calculate_min_radius_see_sphericalpair


##########################
class CylindricaltestShell(testShell):

    self.dr
    self.dz_right
    self.dz_left
    self.scaling_angle_right = None     # this needs to be calculated once

    def determine_possible_shell(self):
        neighbors = self.get_neighbors()
        min_dr_dzright_dzleft = self.get_min_dr_dzright_dzleft()
        max_dr_dzright_dzleft = self.get_max_dr_dzright_dzleft()

        dr_dzright_dzleft = max_dr_dzright_dzleft
        for neighbor in neighbors:

            shell_list = self.get_neighbor_shell_list(neighbor)
            for _, shell_appearance in shell_list:
                new_dr_dzright_dzleft = neighbor.get_dr_dzright_dzleft(shell_appearance, self)

                # if the newly calculated dimensions are smaller than the current one, use them
                dr_dzright_dzleft = min(dr_dzright_dzleft, new_dr_dzright_dzleft)
                if dr_dzright_dzleft < min_dr_dzright_dzleft:
                    dr_dzright_dzleft = None
                    break

        if dr_dzright_dzleft:
            return self.create_new_shell(dr_dzright_dzleft)
        else
            return None

    def get_searchradius():
        return max_that_fits_correctly
            
#####
class PlanarSurfaceSingletestShell(CylindricaltestShell, NonInteractionSingles):

    def get_orientation_vector(self):
        return structure.shape.unit_z   # just copy from structure

    def get_searchpoint(self):
        return particle_pos

    def get_referencepoint(self):
        return particle_pos

    def get_min_dr_dzright_dzleft(self):
        dr = particle_radius * MULTI_SHELL_FACTOR
        dz_right = particle_radius
        dz_left = particle_radius
        return (dr, dz_right, dz_left)

    def get_max_dr_dzright_dzleft(self):
        dr = math.sqrt(geometrycontainer_max**2 - particle_radius**2)
        dz_right = particle_radius
        dz_left = particle_radius
        return (dr, dz_right, dz_left)

    def get_right_scalingcenter(self, orientation_z, orientation_r):
        # returns the scaling center in the coordinate system of the plane through
        # scaling cylinder and the neighboring shell
        return particle_radius * orientation_z

    def get_left_scalingcenter(self, orientation_z, orientation_r):
        return particle_radius * orientation_z

    def get_right_scalingangle(self):
        return Pi/2.0

    def get_left_scalingangle(self):
        return Pi/2.0

#####
class PlanarSurfacePairtestShell(hasCylindricalShell, Others):

    def get_orientation_vector(self):
        return structure.shape.unit_z

    def get_searchpoint(self):
        return CoM

    def get_referencepoint(self):
        return CoM

    def get_min_dr_dzright_dzleft(self):
        do calculation
        return (dr, dz_right, dz_left)

    def get_max_dr_dzright_dzleft(self):
        dr = math.sqrt(geometrycontainer_max**2 - particle_radius**2)
        dz_right = particle_radius
        dz_left = particle_radius
        return (dr, dz_right, dz_left)
        
#####
class CylindricalSurfaceSingletestShell(CylindricaltestShell, NonInteractionSingles):

    def orientation_vector(self):
        return structure.shape.unit_z   # just copy from structure

    def get_searchpoint(self):
        return particle_pos

    def get_referencepoint(self):
        return particle_pos

    def get_min_dr_dzright_dzleft(self):
        dr = particle_radius
        dz_right = particle_radius * MULTI_SHELL_FACTOR
        dz_left = particle_radius * MULTI_SHELL_FACTOR
        return (dr, dz_right, dz_left)

    def get_max_dr_dzright_dzleft(self):
        dr = particle_radius
        dz_right = math.sqrt(geometrycontainer_max**2 - particle_radius**2)
        dz_left = dz_right
        return (dr, dz_right, dz_left)

    def get_right_scalingcenter(self, orientation_z, orientation_r):
        # returns the scaling center in the coordinate system of the plane through
        # scaling cylinder and the neighboring shell
        #
        # Note that the scalingcenter depends on the position of the shell,
        # and therefor orientation_r changes for every shell -> recalculate
        return orientation_r * self.particle_radius

    def get_left_scalingcenter(self, orientation_z, orientation_r):
        # same here
        return orientation_r * self.particle_radius

    def get_right_scalingangle(self):
        return 0

    def get_left_scalingangle(self):
        return 0

#####
class CylindricalSurfacePairtestShell(CylindricaltestShell, Others):

#####
class PlanarSurfaceInteractiontestShell(CylindricaltestShell, Others):

    def get_orientation_vector(self):
        result = do_calculation
        return result

    def get_searchpoint(self):
        result = calculate_searchpoint
        return result

    def get_referencepoint(self):
        return projected_point_of_particle

    def get_min_dr_dzright_dzleft(self):
        dz_right, r = calc_based_on_scaling_stuff
        dz_left = particle_radius
        return (dr, dz_right, dz_left)

    def get_max_dr_dzright_dzleft(self):
        result = calc_based_on_scaling_stuff
        return (dr, dz_right, dz_left)

    def get_right_scalingcenter(self, orientation_z, orientation_r):
        # note calculate this only once since 'orientation_z' doens't change
        # for different shells -> scalingcenter is a fixed point
        h0_right = distance_particle_surface
        return orientation_z * h0_right

    def get_left_scalingcenter(self, orientation_z, orientation_r):
        # note calculate this only once since orientation_z doesn't change
        # for difference shells -> scalingcenter is a fixed point
        return orientation_z * particle_radius

    def get_right_scalingangle(self):
        return bla

    def get_left_scalingangle(self):
        return Pi/2.0

#####
class CylindricalSurfaceInteractiontestShell(CylindricaltestShell, Others):

    def get_orientation_vector(self):
        return structure.shape.unit_z

    def get_searchpoint(self):
        return projected_point_of_particle

    def get_referencepoint(self):
        return projected_point_of_particle

    def get_min_dr_dzright_dzleft(self):
        result = calc_based_on_scaling_stuff
        return (dr, dz_right, dz_left)
        
    def get_max_dr_dzright_dzleft(self):
        result = calc_based_on_scaling_stuff
        return (dr, dz_right, dz_left)

    def get_right_scalingcenter(self, orientation_z, orientation_r):
        # note calculate this only once since 'orientation_z' doens't change
        # for different shells -> scalingcenter is a fixed point
        h0_right = smt
        return orientation_z * h0_right

    def get_left_scalingcenter(self, orientation_z, orientation_r):
        # note calculate this only once since orientation_z doesn't change
        # for difference shells -> scalingcenter is a fixed point
        h0_left = h0_right
        return orientation_z * h0_left

    def get_right_scalingangle(self):
        return bla

    def get_left_scalingangle(self):
        return Pi/2.0

#####
class MixedPair3D2DtestShell(CylindricaltestShell, Others):

    def get_orientation_vector(self):
        result = do_calculation
        return result

    def get_searchpoint(self):
        return caculated_search_point

    def get_referencepoint(self):
        return CoM

    def get_min_dr_dzright_dzleft(self):
        dz_left = particle_radius
        dr, dz_right = calc_based_on_scaling_stuff
        return (dr, dz_right, dz_left)
        
    def get_max_dr_dzright_dzleft(self):
        dz_left = particle_radius
        dr, dz_right = calc_based_on_scaling_stuff
        return (dr, dz_right, dz_left)

    def get_right_scalingcenter(self, orientation_z, orientation_r):
        # note calculate this only once since 'orientation_z' doens't change
        # for different shells
        h0_right = bla
        return orientation_z * h0_right

    def get_left_scalingcenter(self, orientation_z, orientation_r):
        # note calculate this only once since orientation_z doesn't change
        # for difference shells
        return orientation_z * particle_radius

    def get_right_scalingangle(self):
        return bla

    def get_left_scalingangle(self):
        return Pi/2.0

#####
#class MixedPair3D1DtestShell(CylindricaltestShell, Others):
