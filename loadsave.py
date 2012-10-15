#!/usr/env python

import sys
import math
import numpy

from array import *

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
from gfrdbase import DomainEvent
import model
import _gfrd

__all__ = [ 'save_state', 'load_state' ]



def save_state(simulator, filename):
    

    #### BURST ALL DOMAINS ####
    # This is an important precondition for this
    # save routine to work properly.
    # In particular running this without bursting
    # all the domains will mess up the EventScheduler
    # because it will reconstruct it with different
    # event IDs
    simulator.burst_all_domains()

    # Re-seed the random number generator
    # (with a randomly drawn seed)
    # This seed has to be saved and used to
    # (re)seed the newly setup simulator when
    # loading to ensure that we will reproduce
    # the same trajectory.
    new_seed = simulator.rng.uniform_int(0, 1000000)
    simulator.reset_seed(new_seed)

    # Reset counters
    N_structure_types = 0
    N_species         = 0
    N_rules           = 0

    #### DEFINE THE CONFIG PARSER OBJECT ####
    cp = CP.ConfigParser()
    cp.optionxform = str # for upper case option names

    # Now add the relevant info in separate sections
    #### MODEL, WORLD and SEED ####
    cp.add_section('MODEL')
    cp.set('MODEL', 'world_size', simulator.world.world_size)

    cp.add_section('WORLD')
    cp.set('WORLD', 'world_size', simulator.world.world_size)
    cp.set('WORLD', 'matrix_size', simulator.world.matrix_size)

    cp.add_section('SEED')
    cp.set('SEED', 'seed', new_seed)

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

        # Update counter
        N_structure_types += 1

    #### SPECIES ####
    for species in simulator.get_species():

        id_int = id_to_int(species.id)
        sectionname = 'SPECIES_' + str(id_int)
        cp.add_section(sectionname)

        cp.set(sectionname, 'id', id_int)
        #cp.set(sectionname, 'name', species.name]) # TODO SpeciesInfo does not have this property for some reason...
        cp.set(sectionname, 'radius', species.radius)
        cp.set(sectionname, 'D', species.D)
        cp.set(sectionname, 'v', species.v)
        cp.set(sectionname, 'structure_type_id', id_to_int(species.structure_type_id))

        # Update counter
        N_species += 1

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
        for i, name in structure_type_names: ## FIXME replace by dictionary?
            if i == sid_int:
                st_name = name

        # Find the keyword for the geometric type of this structure
        for object_type, key in structure_keywords:
            if isinstance(structure, object_type):
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

    #### PARTICLE ORDER IN SCHEDULER ####
    # Pop the events from the scheduler and remember
    # them to push them back again afterwards.
    # Event IDs will change but this should not matter
    # for the simulated trajectory as such.    
    eventlist = []
    scheduler_order = []
    while not simulator.scheduler.size == 0 :

        id, event = simulator.scheduler.pop()
        domain = simulator.domains[event.data]

        # Store the events in a list
        # inserting from the beginning; thus it
        # will already have the right order for
        # reinsertion
        eventlist.insert(0, (event, domain))

        pid, particle = domain.pid_particle_pair        
        pid_int = id_to_int(pid)

        scheduler_order.append(pid_int)
    
    # Just to be sure...
    assert simulator.scheduler.size == 0
    # Put the events back into the scheduler in reverse order
    for (event, domain) in list(eventlist):
        event_id = simulator.scheduler.add(DomainEvent(event.time, domain))

    # Create a new section in the save file
    sectionname = 'SCHEDULER'
    cp.add_section(sectionname)
    cp.set(sectionname, 'particle_order', list(scheduler_order))    

    #### ADD COUNTERS TO MODEL SECTION ####
    cp.set('MODEL', 'N_structure_types', N_structure_types)
    cp.set('MODEL', 'N_species', N_species)

    #### WRITE FILE ####
    with open(filename, 'wb') as outfile:
        cp.write(outfile)



def load_state(filename):
    """ Loads a file previously generated with save_state() and
        constructs a model, world and info needed to construct the 
        EGFRDSimulator.       
    """

    #### DEFINE THE CONFIG PARSER OBJECT ####
    cp = CP.ConfigParser()
    cp.optionxform = str # for upper case option names

    # Load the saved file
    # TODO Check for right format
    cp.read(filename)

    print 'Sections in file : ' + str(cp.sections()) ### TESTING

    #### FIRST GET GLOBAL INFO ####
    world_size  = cp.getfloat('WORLD', 'world_size')
    matrix_size = cp.getint('WORLD', 'matrix_size')
    N_structure_types = cp.getint('MODEL', 'N_structure_types')
    N_species   = cp.getint('MODEL', 'N_species')
    seed        = cp.getint('SEED', 'seed')

    #### CREATE THE WORLD AND THE MODEL ####
    m = model.ParticleModel(world_size)


    #### STRUCTURE_TYPES ####
    structure_types_dict = {}  # will map the old (read-in) ID to the StructureType

    structure_types_sections = filter_sections(cp.sections(), 'STRUCTURETYPE')
    # Now read in the the structure types and add them to the new model
    # in the right order.
    # Note that add_structure_type() will create IDs automatically,
    # so it is important that we sort the sectionnames list by the ids.
    for sectionname in sorted(structure_types_sections, key = lambda name : name_to_int(name)):
        
        id      = cp.getint(sectionname, 'id')
        name    = cp.get(sectionname, 'name')        

        # Assert that the read id corresponds to the one in the sectionname
        assert name_to_int(sectionname) == id

        # Create a new StructureType object and add it to the model
        structure_type = _gfrd.StructureType()
        structure_type['name'] = name

        if id > 1:  # Skip adding the default structure_type a second time
            m.add_structure_type(structure_type)
            assert id_to_int(structure_type.id) == id
            structure_types_dict[id] = structure_type
            print 'Added ' + str(structure_type) ### TESTING

    print 'structure_types_dict = ' + str(structure_types_dict) ### TESTING


    #### SPECIES ####
    # Now the same for the species
    # Again we have to make sure we add the species in the right order => sort sections list
    species_dict = {}  # will map the old (read-in) ID to the Species

    species_sections = filter_sections(cp.sections(), 'SPECIES')
    for sectionname in sorted(species_sections, key = lambda name : name_to_int(name)):
        
        id      = cp.getint(sectionname, 'id')
        #name    = cp.getfloat(sectionname, 'name') # FIXME
        radius  = cp.getfloat(sectionname, 'radius')
        D       = cp.getfloat(sectionname, 'D')
        v       = cp.getfloat(sectionname, 'v')
        structure_type_id = cp.getint(sectionname, 'structure_type_id')

        # Assert that the read id corresponds to the one in the sectionname
        assert name_to_int(sectionname) == id

        # FIXME name cannot be output yet, so we have to make up one here
        name = 'Species_'+str(id)

        # Get the structure type for this species from the structure_type_id
        if structure_type_id > 1:
            structure_type = structure_types_dict[structure_type_id]
            assert id_to_int(structure_type.id) == structure_type_id

        # Create a new Species object and add it to the model
        if structure_type_id == 1:
            species = model.Species(name, D, radius)
        elif v == 0:
            species = model.Species(name, D, radius, structure_type)
        else:
            species = model.Species(name, D, radius, structure_type, v)

        m.add_species_type(species)
        assert id_to_int(species.id) == id
        species_dict[id] = species
        print 'Added ' + str(species) ### TESTING

    print 'species_dict = ' + str(species_dict) ### TESTING


    #### RULES ####
    rules_dict = {}  # will map the old (read-in) ID to the reaction rule

    rules_sections = filter_sections(cp.sections(), 'REACTIONRULE')
    for sectionname in sorted(rules_sections, key = lambda name : name_to_int(name)):

        rtype        = cp.get(sectionname, 'type')
        rate         = cp.getfloat(sectionname, 'rate')
        reactant_ids = eval(cp.get(sectionname, 'reactant_ids')) # eval makes sure we get a list
        product_ids  = eval(cp.get(sectionname, 'product_ids'))

        rule = None
        if rtype == 'unimolecular':

            assert len(reactant_ids) == 1 and reactant_ids[0] != None and \
                   len(product_ids)  == 1 and product_ids[0]  != None, \
                     'Saved reaction rule does not have the proper signature: rule = %s' % sectionname

            # Get the involved species
            reactant = species_dict[int(reactant_ids[0])]
            product  = species_dict[int(product_ids[0])]

            # Create the rule
            rule = model.create_unimolecular_reaction_rule(reactant, product, rate)

        elif rtype == 'unbinding':

            assert len(reactant_ids) == 1 and reactant_ids[0] != None and \
                   product_ids[0] != None and product_ids[1]  != None, \
                     'Saved reaction rule does not have the proper signature: rule = %s' % sectionname

            # Get the involved species
            reactant = species_dict[int(reactant_ids[0])]
            product0 = species_dict[int(product_ids[0])]
            product1 = species_dict[int(product_ids[1])]

            # Create the rule
            rule = model.create_unbinding_reaction_rule(reactant, product0, product1, rate)

        elif rtype == 'binding_particle':
            
            assert reactant_ids[0] != None and reactant_ids[1] != None and \
                   len(product_ids) == 1   and product_ids[0]  != None, \
                     'Saved reaction rule does not have the proper signature: rule = %s' % sectionname

            # Get the involved species
            reactant0 = species_dict[int(reactant_ids[0])]
            reactant1 = species_dict[int(reactant_ids[1])]
            product   = species_dict[int(product_ids[0])]

            # Create the rule
            rule = model.create_binding_reaction_rule(reactant0, reactant1, product, rate)

        elif rtype == 'binding_surface':
            
            assert reactant_ids[0] != None and reactant_ids[1] != None and \
                   len(product_ids) == 1   and product_ids[0]  != None, \
                     'Saved reaction rule does not have the proper signature: rule = %s' % sectionname

            # Get the involved species
            structure_type = structure_types_dict[int(reactant_ids[0])]
            reactant       = species_dict[int(reactant_ids[1])]
            product        = species_dict[int(product_ids[0])]

            # Create the rule
            rule = model.create_binding_reaction_rule(reactant, structure_type, product, rate)

        else:

            raise LoadSaveError('Cannot create reaction rule when loading state: unknown reaction rule type.')

        # If we have successfully constructed a reaction rule, add it to the model
        # and to the internal rules dictionary
        if rule:
            
            m.add_reaction_rule(rule)
            rules_dict[name_to_int(sectionname)] = rule
            print 'Added ' + str(rule) + ', type = %s' % str(rtype) ### TESTING

    print 'rules_dict = ' + str(rules_dict) ### TESTING

    #### CREATE THE WORLD ####
    w = create_world(m, matrix_size)
    # This will come in handy later
    def_structure_id = w.get_def_structure_id()

    #### STRUCTURES ####
    structures_dict = {}  # will map the old (read-in) ID to the structure

    structures_sections = filter_sections(cp.sections(), 'STRUCTURE')
    for sectionname in sorted(structures_sections, key = lambda name : name_to_int(name)):

        # Get the geometry key
        # TODO This assumes that the section name has the right signature
        #      => Implement a check
        geometry_key = sectionname.split(get_default_separator())[-2]

        # Find the geometric type of this structure according to its key
        for object_type, key in structure_keywords:
            if geometry_key == key:
                structure_object_type = object_type

        # Read in the information common to all object types
        id        = cp.getint(sectionname, 'id')        
        name      = cp.get(sectionname, 'name')
        parent_id = cp.getint(sectionname, 'parent_id')
        st_id     = cp.getint(sectionname, 'structure_type_id')        
        st_name   = cp.get(sectionname, 'structure_type_name')
        position  = eval(cp.get(sectionname, 'position'))
                    # eval makes sure we import a list

        # Assert that the read id corresponds to the one in the sectionname
        # This is particularly important in this part of the loading routine
        # because structures may be defined as substructures of previously
        # defined structures. Thus we want to define them in the order in
        # which they were defined before they have been saved to the file.
        assert name_to_int(sectionname) == id

        parent_structure = w.get_structure(def_structure_id) # TESTING

        structure = None
        if structure_object_type == CuboidalRegion:
        
            unit_x      = eval(cp.get(sectionname, 'unit_x'))
            unit_y      = eval(cp.get(sectionname, 'unit_y'))
            unit_z      = eval(cp.get(sectionname, 'unit_z'))            
            half_extent = eval(cp.get(sectionname, 'half_extent'))

            assert id == id_to_int(def_structure_id)
            # TODO Add more than default cuboidal regions
            # Right now we assume that there is only one which
            # is created automatically via create_world()            

        elif structure_object_type == SphericalSurface:
            
            radius = cp.getfloat(sectionname, 'radius')

            # Create the structure
            structure = model.create_spherical_surface(structure_type.id, name, position, radius, parent_structure.id)

        elif structure_object_type == CylindricalSurface:
            
            unit_z      = eval(cp.get(sectionname, 'unit_z'))            
            half_length = cp.getfloat(sectionname, 'half_length')
            radius      = cp.getfloat(sectionname, 'radius')

        elif structure_object_type == DiskSurface:
            
            unit_z      = eval(cp.get(sectionname, 'unit_z'))            
            radius      = cp.getfloat(sectionname, 'radius')

            structure_type = structure_types_dict[st_id]            

            # Create the structure
            structure = model.create_disk_surface(structure_type.id, name, position, radius, unit_z, parent_structure.id)

        elif structure_object_type == PlanarSurface:

            unit_x      = eval(cp.get(sectionname, 'unit_x'))
            unit_y      = eval(cp.get(sectionname, 'unit_y'))
            unit_z      = eval(cp.get(sectionname, 'unit_z'))            
            half_extent = eval(cp.get(sectionname, 'half_extent'))

        # Add it to the world
        if structure:
            w.add_structure(structure)
            #assert id_to_int(structure.id) == id # TODO
            structures_dict[id] = structure
            print 'Added ' + str(structure) + ', type = %s, parent = %s' % (structure_type, parent_id) ### TESTING

    print 'structures_dict = ' + str(structures_dict) ### TESTING


##########################
#### HELPER FUNCTIONS ####
##########################
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


def get_default_separator():

    return '_'


def separate(sectionname):

    separator = get_default_separator()
    part = sectionname.partition(separator)

    return (part[0], part[2])

    # TODO Check for right tag format
    # raise LoadSaveError('Section name in input file did not have the required format.')


def filter_sections(sectionlist, tagstring):

    assert isinstance(tagstring, str)

    return [section for section in sectionlist if separate(section)[0] == tagstring]


def name_to_int(sectionname):

    separator = get_default_separator()

    return int(sectionname.split(separator)[-1])


structure_keywords = [
        (CuboidalRegion,     'CUBE'),
        (SphericalSurface,   'SPHERE'),
        (CylindricalSurface, 'CYLINDER'),
        (DiskSurface,        'DISK'),
        (PlanarSurface,      'PLANE')
        ]


class LoadSaveError(Exception):
    pass

