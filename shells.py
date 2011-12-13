#!/usr/bin/python

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
    def get_dr_dzright_dzleft_to_shell(self, shell_appearance, domain_to_scale):
        # This will scale the 'domain_to_scale' (a cylinder) using the 'shell_appearance' as the limiting shell
        # CylindricaltestShell ('domain_to_scale') -> SphericalShell ('self' with 'shell_appearance')

        assert type(shell_appearance, Sphere)        # because the shell actually originated here
        # Do Laurens' algorithm (part1)

        return (r, z_right, z_left)


    def get_radius_to_shell(self, shell_appearance, domain_to_scale):
        # SphericaltestShell ('domain_to_scale') -> SphericalShell ('self' with 'shell_appearance')

        assert type(shell_appearance, Sphere)        # because the shell actually originated here
        # Do simple distance calculation to sphere

        return r

#####
class SphericalSingle(hasSphericalShell):

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

# Others
#####
class SphericalPair(hasSphericalShell):

    def shell_list_for_single(self):
        return self.shell_list
    def shell_list_for_other(self):
        return self.shell_list

#####
class Multi(hasSphericalShell):

    def shell_list_for_single(self):
        return self.shell_list
    def shell_list_for_other(self):
        return self.shell_list



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
class PlanarSurfaceSingle(hasCylindricalShell):

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
class CylindricalSurfaceSingle(hasCylindricalShell):

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

# Others
#####
class PlanarSurfacePair(hasCylindricalShell):

    def shell_list_for_single(self):
        return self.shell_list
    def shell_list_for_other(self):
        return self.shell_list

#####
class PlanarSurfaceInteraction(hasCylindricalShell):

    # these return uncorrected dimensions because Interactions seem the same
    def shell_list_for_single(self):
        return self.shell_list
    def shell_list_for_other(self):
        return self.shell_list

#####
class CylindricalSurfacePair(hasCylindricalShell):
class CylindricalSurfaceInteraction(hasCylindricalShell):
class MixedPair3D2D(hasCylindricalShell):
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
        return shell_container.get_max()

#####
class SphericalSingletestShell(testSphericalShell):

    self.center = particle_position

    def get_min_radius(self):
        return particle_radius * MULTI_SHELL_FACTOR     # the minimum radius of a NonInteractionSingle

    # always the same for singles
    def get_neighbor_shell_list(self, neighbor):
        # note that this returns the list of shells, NOT dr, dz_right, dz_left
        return neighbor.shell_list_for_single()

#####
class SphericalPairtestShell(testSphericalShell):

    self.center = calculate CoM
    def get_min_radius:
        return calculate_min_radius_see_sphericalpair

    # always the same for other (InteractionSingle, Pair, Multi)
    def get_neighbor_shell_list(self, neighbor):
        # note that this returns the list of shells, NOT dr, dz_right, dz_left
        return neighbor.shell_list_for_other()


##########################
class CylindricaltestShell(testShell):

    self.dr
    self.dz_right
    self.dz_left

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
class PlanarSurfaceSingletestShell(CylindricaltestShell):

    def get_orientation_vector(self):
        return structure.shape.unit_z   # just copy from structure

    def get_searchpoint(self):
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

    # always the same for NonInteractionSingle
    def get_neighbor_shell_list(self, neighbor):
        # note that this returns the list of shells, NOT dr, dz_right, dz_left
        return neighbor.shell_list_for_single()

#####
class PlanarSurfaceInteractiontestShell(CylindricaltestShell):

    def get_orientation_vector(self):
        result = do_calculation
        return result

    def get_searchpoint(self):
        result = calculate_searchpoint
        return result

    def get_min_dr_dzright_dzleft(self):
        dz_right, r = calc_based_on_scaling_stuff
        dz_left = particle_radius
        return (dr, dz_right, dz_left)

    def get_max_dr_dzright_dzleft(self):
        result = calc_based_on_scaling_stuff
        return (dr, dz_right, dz_left)

    # always the same for other (InteractionSingle, Pair, Multi)
    def get_neighbor_shell_list(self, neighbor):
        # note that this returns the list of shells, NOT dr, dz_right, dz_left
        return neighbor.shell_list_for_other()

