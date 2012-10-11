#!/usr/env python

import math
import numpy

import ConfigParser as CP

from _gfrd import (
    Event,
    EventScheduler,
    Particle,
    SphericalShell,
    CylindricalShell,
    DomainIDGenerator,
    ShellIDGenerator,
    DomainID,
    ParticleContainer,
    CuboidalRegion,
    SphericalSurface,
    CylindricalSurface,
    DiskSurface,
    PlanarSurface,
    Surface,
    _random_vector,
    Cylinder,
    Disk,
    Sphere,
    Plane,
    NetworkRulesWrapper,
    )

from gfrdbase import *

__all__ = [ 'save_state', 'load_state' ]



def save_state(simulator, filename):
    

    #### DEFINE THE CONFIG PARSER OBJECT ####
    cp = CP.ConfigParser()
    cp.optionxform = str # for upper case option names

    # Now add the relevant info in separate sections
    #### MODEL and WORLD ####
    cp.add_section('MODEL')
    cp.set('MODEL', 'world_size', simulator.world.world_size)

    cp.add_section('WORLD')
    cp.set('WORLD', 'world_size', simulator.world.world_size)
    cp.set('WORLD', 'matrix_size', simulator.world.matrix_size)

    #### SPECIES ####
    for species in simulator.get_species():

        id_int = id_to_int(species.id)
        sectionname = 'SPECIES_' + str(id_int)
        cp.add_section(sectionname)

        cp.set(sectionname, 'id', id_int)
        cp.set(sectionname, 'radius', species.radius)
        cp.set(sectionname, 'D', species.D)
        cp.set(sectionname, 'v', species.v)
        cp.set(sectionname, 'structure_type_id', id_to_int(species.structure_type_id))

    #### STRUCTURE TYPES ####
    structure_type_names = [] # will store the names in a lookup list
    for structure_type in simulator.get_structure_types():

        id_int = id_to_int(structure_type.id)

        if structure_type['name'] != '' :
            name = structure_type['name']

        elif id_int == 1 : # the std. structure type
            name = 'default'

        sectionname = 'STRUCTURETYPE_' + str(id_int)
        cp.add_section(sectionname)        

        cp.set(sectionname, 'id', id_int)
        cp.set(sectionname, 'name', name)

        # Remember name for later
        structure_type_names.append( (id_int, name) )

    #### RULES ####
    # TODO ...
    for species in simulator.get_species():

        rules = simulator.network_rules.query_reaction_rule(species)
        

    #### STRUCTURES ####
    for structure in simulator.get_structures():

        id_int  = id_to_int(structure.id)
        sid_int = id_to_int(structure.sid)
        parent_id_int = id_to_int(structure.structure_id)
            # the structure_id property of structures is the parent structure id

        # Find the name of the structure type
        for i, name in structure_type_names:
            if i == sid_int:
                st_name = name

        # Find the keyword for the geometric type of this structure
        for instance, key in structure_keywords:
            if isinstance(structure, instance):
                geometry_key = key

        sectionname = 'STRUCTURE_' + geometry_key + '_' + str(id_int)
        cp.add_section(sectionname)

        cp.set(sectionname, 'id', id_int)
        cp.set(sectionname, 'parent_id', parent_id_int)
        cp.set(sectionname, 'structure_type_id', sid_int)
        cp.set(sectionname, 'structure_type_name', st_name)
        cp.set(sectionname, 'name', structure.name)

        cp.set(sectionname, 'position', list(structure.shape.position))            

        if isinstance(structure, CuboidalRegion):
            
            cp.set(sectionname, 'unit_x', list(structure.shape.unit_x))
            cp.set(sectionname, 'unit_y', list(structure.shape.unit_y))
            cp.set(sectionname, 'unit_z', list(structure.shape.unit_z))
            cp.set(sectionname, 'half_extent', list(structure.shape.half_extent))

        elif isinstance(structure, SphericalSurface):

            cp.set(sectionname, 'radius', structure.shape.radius)

        elif isinstance(structure, CylindricalSurface):

            cp.set(sectionname, 'unit_z', list(structure.shape.unit_z))
            cp.set(sectionname, 'radius', structure.shape.radius)
            cp.set(sectionname, 'half_length', structure.shape.half_length)

        elif isinstance(structure, DiskSurface):

            cp.set(sectionname, 'unit_z', list(structure.shape.unit_z))
            cp.set(sectionname, 'radius', structure.shape.radius)

        elif isinstance(structure, PlanarSurface):

            cp.set(sectionname, 'unit_x', list(structure.shape.unit_x))
            cp.set(sectionname, 'unit_y', list(structure.shape.unit_y))
            cp.set(sectionname, 'unit_z', list(structure.shape.unit_z))
            cp.set(sectionname, 'half_extent', list(structure.shape.half_extent))


    #### STRUCTURE LINKS ####
    # TODO

    #### PARTICLES ####
    pid_particle_pairs = list(simulator.world)
    pid_particle_pairs.sort()

    for pp in pid_particle_pairs:

        pid = pp[0]
        particle = pp[1]

        species_id_int   = id_to_int(particle.sid)
        structure_id_int = id_to_int(pp[1].structure_id)

        pid_int = id_to_int(pid)

        sectionname = 'PARTICLE_' + str(pid_int)
        cp.add_section(sectionname)

        cp.set(sectionname, 'id', pid_int)
        cp.set(sectionname, 'species_id', species_id_int)
        cp.set(sectionname, 'radius', particle.radius)
        cp.set(sectionname, 'D', particle.D)
        cp.set(sectionname, 'v', particle.v)
        cp.set(sectionname, 'position', list(particle.position))
        cp.set(sectionname, 'structure_id', structure_id_int)


    #### WRITE FILE ####
    with open(filename, 'wb') as outfile:
        cp.write(outfile)



def load_state(simulator, filename):

    pass


def id_to_int(ID):
    """ Strips an eGFRD-type ID-specifier off everything
        but the integer ID numbers.

        I.e., PID(0:2) for example will be converted to 2.

        TODO Is the first number entry not important, too?
    """

    p = str(ID).partition(':')

    if p[2] != '':
        q = p[2].partition(')')

    if q[0] != '':
        return int(q[0])

    else:
        raise RuntimeError('Could not extract number, probably the argument is not a valid ID.')


structure_keywords = [
        (CuboidalRegion,     'CUBE'),
        (SphericalSurface,   'SPHERE'),
        (CylindricalSurface, 'CYLINDER'),
        (DiskSurface,        'DISK'),
        (PlanarSurface,      'PLANE')
        ]