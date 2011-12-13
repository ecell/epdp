#!/usr/bin/python

class hasShell(object):

class hasSphericalShell(hasShell):

    def __init__(self, sphericaltestShell):
        self.center = sphericaltestShell.center
        self.radius = sphericaltestShell.radius

    # these return potentially corrected dimensions
    def get_center:
    def get_radius:

class SphericalSingle(hasSphericalShell):

    def get_radius_for_single(self):
        return max(self.get_min_radius(), self.radius)
    def get_radius_for_other(self):
        return max(particle_radius * SINGLE_SHELL_FACTOR, self.radius)

class SphericalPair(hasSphericalShell):

    def get_radius_for_single(self):
        return self.radius
    def get_radius_for_other(self):
        return self.radius

class Multi(hasSphericalShell):




class hasCylindricalShell(hasShell):

    def __init__(self, cylindricaltestShell):
        self.center = uitrekenen zie interactions
        self.radius = cylindricaltestShell.dr
        self.half_length = uitrekenen zie interactions

    def get_center:

class PlanarSurfaceSingle(hasCylindricalShell):

    # these return corrected dimensions
    def get_dimensions_for_single(self):
        return max(self.get_min_dimensions, (self.radius, self.half_length))
    def get_dimensions_for_other(self):
        return (max(particle_radius * SINGLE_SHELL_FACTOR, self.radius), self.half_length)

class PlanarSurfaceInteraction(hasCylindricalShell):

    # these return corrected dimensions
    def get_dimensions_for_single(self):
        return (self.radius, self.half_length)
    def get_dimensions_for_other(self):
        return (self.radius, self.half_length)




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
        for neighbor in neighbors:
            if isinstance(neighbor, hasCylindricalShell):
                dr_dzright_dzleft = self.get_dr_dzright_dzleft_CylindricalNeighbor(neighbor)
            elif isinstance(neighbor, hasSphericalShell):
                dr_dzright_dzleft = self.get_dr_dzright_dzleft_SphericalNeighbor(neighbor)
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
            
    def get_dr_dzright_dzleft_CylindricalNeighbor(self, neighbor):
        # get's the maximal (dr, dz_right, dz_left) values for a neighboring CYLINDRICAL domain
        # CylindricaltestShell (self) -> CylindricalShell (neighbor)
        # Laurens' algorithm (part2)
        neighbor_dimensions = self.get_neighbor_dimensions(neighbor)
        return (r, z_right, z_left)

    def get_dr_dzright_dzleft_SphericalNeighbor(self, neighbor):
        # get's the maximal (dr, dz_right, dz_left) values for a neighboring SPHERICAL domain
        # CylindricaltestShell (self) -> SphericalShell (neighbor)
        # Laurens' algorithm (part1)
        neighbor_radius = self.get_neighbor_radius(neighbor)
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
    def get_neighbor_dimensions(self, neighbor):
        # note that this returns (radius, half_length) tuple, NOT dr, dz_right, dz_left
        return neighbor.get_dimensions_for_single()
    def get_neighbor_radius(self, neighbor):
        return neighbor.get_radius_for_single()


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
    def get_neighbor_dimensions(self, neighbor):
        # note that this returns (radius, half_length) tuple, NOT dr, dz_right, dz_left
        return neighbor.get_dimensions_for_other()
    def get_neighbor_radius(self, neighbor):
        return neighbor.get_radius_for_other()

