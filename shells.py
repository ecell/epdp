#!/usr/bin/python

class hasShell(object):

    def shell_list(self):
            return self.shell_list

class hasSphericalShell(hasShell):

    self.shell

    def __init__(self, sphericaltestShell):
        self.center = sphericaltestShell.center
        self.radius = sphericaltestShell.radius

    # these return potentially corrected dimensions
    def get_center:
    def get_radius:

class SphericalSingle(hasSphericalShell):

    def shell_list_for_single(self):
        if self.get_min_radius() > self.shell.radius:
            fake_shell = make_new_shell_with_min_radius
            return [fake_shell, ]
        else:
            return real_shell_list

    def shell_list_for_other(self):
        if shellradius < particle_radius * SINGLE_SHELL_FACTOR:
            fake_shell = make_new_shell (particle_radius * SINGLE_SHELL_FACTOR)
            return [fake_shell, ]
        else:
            return real_shell_list


class SphericalPair(hasSphericalShell):

    def shell_list_for_single(self):
        return real_shell_list
    def shell_list_for_other(self):
        return real_shell_list

class Multi(hasSphericalShell):

    def shell_list_for_single(self):
        return self.shell_list
    def shell_list_for_other(self):
        return self.shell_list


class hasCylindricalShell(hasShell):

    def __init__(self, cylindricaltestShell):
        self.center = uitrekenen zie interactions
        self.radius = cylindricaltestShell.dr
        self.half_length = uitrekenen zie interactions

    def get_center:

class PlanarSurfaceSingle(hasCylindricalShell):

    # these return corrected dimensions, since we reserve more space for the NonInteractionSingle
    def shell_list_for_single(self):
        if self.get_min_dimensions > shell.dimensions:
            fake_shell = create_shell(min_dimensions)
            return [fake_shell, ]
        else:
            return self.shell_list

    def shell_list_for_other(self):
        return (max(particle_radius * SINGLE_SHELL_FACTOR, self.radius), self.half_length)

class PlanarSurfaceInteraction(hasCylindricalShell):

    # these return uncorrected dimensions because Interactions seem the same
    def shell_list_for_single(self):
        return self.shell_list
    def shell_list_for_other(self):
        return self.shell_list




#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
class testShell(object):

    def get_neighbors(self):
        searchpoint = self.get_searchpoint()
        radius = self.get_searchradius()
        return self.geometrycontainer.get_neigbors(searchpoint, radius)

##########################
class SphericaltestShell(testShell):

    self.radius
    self.center

    def create_shell:
        neighbors = self.get_neighbors()
        min_radius = self.get_min_radius()
        max_radius = self.get_max_radius()
        for neighbor in neighbors:
            if isinstance(neighbor, hasCylindricalShell):
                neighbor_dimensions = self.get_neighbor_dimensions(neighbor)
                radius = self.scale_max_sphere_cylinder(neighbor_dimensions)
            elif isinstance(neighbor, hasSphericalShell):
                neighbor_radius = self.get_neighbor_radius(neighbor)
                radius = self.scale_max_sphere_sphere(neighbor_radius)
            min(radius, max_radius)
            if radius < min_radius:
                break
        create_new_shell(radius)
            

    def get_searchpoint(self):
        return self.center
    def get_searchradius(self):
        return self.geometrycontainer.max()

    def get_max_radius(self):
        return shell_container.get_max()

    def get_radius_SphericalNeighbor(self, neighbor_radius):
        # SphericaltestShell (self) -> SphericalShell (neighbor)
        # Simple distance calculation to sphere
        return r

    def get_radius_CylindricalNeighbor(self, neighbor_dimensions):
        # SphericaltestShell (asker) -> CylindricalShell (self)
        # Simple distance calculation to cylinder
        return r

class SphericalSingletestShell(testSphericalShell):

    self.center = particle_position
    def get_min_radius(self):
        return particle_radius * MULTI_SHELL_FACTOR     # the minimum radius of a NonInteractionSingle

class SphericalPairtestShell(testSphericalShell):

    self.center = calculate CoM
    def get_min_radius:
        return calculate_min_radius_see_sphericalpair



###########################
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
            for _, shell in shell_list:
                # use the right algorithm for the neighbor
                if isinstance(shell, CylindricalShell):
                    dr_dzright_dzleft = self.get_dr_dzright_dzleft_CylindricalNeighbor(shell)
                elif isinstance(shell, SphericalShell):
                    dr_dzright_dzleft = self.get_dr_dzright_dzleft_SphericalNeighbor(shell)

            # if the newly calculated dimensions are smaller than the current one, use them
            min(dr_dzright_dzleft, max_dr_dzright_dzleft)
            if dr_dzright_dzleft < min_dr_dzright_dzleft:
                dr_dzright_dzleft = None
                break
        if dr_dzright_dzleft:
            return self.create_new_shell(dr_dzright_dzleft)
        else
            return None

    def get_searchradius():
        return max_that_fits_correctly
            
    def get_dr_dzright_dzleft_CylindricalNeighbor(self, shell):
        # get's the maximal (dr, dz_right, dz_left) values for a neighboring CYLINDRICAL domain
        # CylindricaltestShell (self) -> CylindricalShell (neighbor)
        # Laurens' algorithm (part2)
        assert type(shell, Cylinder)

        return (r, z_right, z_left)

    def get_dr_dzright_dzleft_SphericalNeighbor(self, neighbor):
        # get's the maximal (dr, dz_right, dz_left) values for a neighboring SPHERICAL domain
        # CylindricaltestShell (self) -> SphericalShell (neighbor)
        # Laurens' algorithm (part1)
        assert type(shell, Sphere)

        return (r, z_right, z_left)

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

