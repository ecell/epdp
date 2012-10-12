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
    extracted_rules = []  # to avoid double-extraction
    for species0 in simulator.get_species():

        # Reactions with one reactant (unimolecular or unbinding)
        rules = simulator.network_rules.query_reaction_rule(species0)
        
        for rr in rules:

              if not rr.id in extracted_rules:

                  rate = rr.k

                  r0 = rr.reactants[0]
                  reactant_ids = ( id_to_int(r0), )

                  if len(rr.products) == 1:

                      rtype = 'unimolecular'
                      p0 = rr.products[0]
                      product_ids  = ( id_to_int(p0), )

                  elif len(rr.products) == 2:

                      rtype = 'unbinding'
                      p0 = rr.products[0]
                      p1 = rr.products[1]
                      product_ids  = ( id_to_int(p0), id_to_int(p1) )

                  else:
                      raise LoadSaveError('Unimolecular or unbinding reaction rule has more than two products. Something is utterly wrong.')

                  sectionname = 'REACTIONRULE_' + str(rr.id)
                  cp.add_section(sectionname)

                  cp.set(sectionname, 'type', rtype)
                  cp.set(sectionname, 'reactant_ids', reactant_ids)
                  cp.set(sectionname, 'product_ids', product_ids)
                  cp.set(sectionname, 'rate', rate)

                  extracted_rules.append(rr.id)

        # Reactions with two reactants (particle species)
        for species1 in list(simulator.get_species()) + list(simulator.get_structure_types()):
              
            rules = simulator.network_rules.query_reaction_rule(species0, species1)

            for rr in rules:

                if not rr.id in extracted_rules and not len(rr.products)==0:

                    rate = rr.k

                    # Check whether we have a surface binding reaction.
                    # Since only species1 loops over StructureTypes we only have to check
                    # whether species1 is a SstructureType or not.
                    if species1.id in [st.id for st in list(simulator.get_structure_types())]:
                        # We must construct an id list from the iterator because otherwise
                        # the comparison does not work
                        rtype = 'binding_surface'

                    else:
                        rtype = 'binding_particle'

                    # Make sure that the first reactant stores species1, which is the surface
                    # in case we have a surface reaction.
                    if rr.reactants[0] == species1.id:
                        r0 = rr.reactants[0]
                        r1 = rr.reactants[1]
                    else:
                        # reverse order
                        r0 = rr.reactants[1]
                        r1 = rr.reactants[0]

                    reactant_ids = ( id_to_int(r0), id_to_int(r1) )

                    # Check whether there is a legal number of products
                    if len(rr.products) == 0:

                        product_ids  = None
                        rtype = 'annihilation'

                    elif len(rr.products) == 1:

                        p0 = rr.products[0]
                        product_ids  = ( id_to_int(p0), )

                    else:
                        raise LoadSaveError('Binding reaction rule has %s product(s). Something is utterly wrong.' % str(len(rr.products)) )           

                    # Construct the section info
                    sectionname = 'REACTIONRULE_' + str(rr.id)
                    cp.add_section(sectionname)

                    cp.set(sectionname, 'type', rtype)
                    cp.set(sectionname, 'reactant_ids', reactant_ids)
                    cp.set(sectionname, 'product_ids', product_ids)
                    cp.set(sectionname, 'rate', rate)

                    extracted_rules.append(rr.id)

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


    #### STRUCTURE CONNECTIVITY ####
    for structure in simulator.get_structures():

        if isinstance(structure, PlanarSurface):

            id_int  = id_to_int(structure.id)

            sectionname = 'STRUCTURECONNECTION_' + str(id_int)
            cp.add_section(sectionname)

            # Retrieve the neighbor structure IDs for this structure
            for n in range(0, 4):

                nid = simulator.world.get_neighbor_id(structure, n)
                nid_int  = id_to_int(nid)

                cp.set(sectionname, 'neighbor_id_'+str(n), nid_int)

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
        raise LoadSaveError('Could not extract number, probably the argument is not a valid ID.')


structure_keywords = [
        (CuboidalRegion,     'CUBE'),
        (SphericalSurface,   'SPHERE'),
        (CylindricalSurface, 'CYLINDER'),
        (DiskSurface,        'DISK'),
        (PlanarSurface,      'PLANE')
        ]


class LoadSaveError(Exception):
    pass

