#ifndef WORLD_HPP
#define WORLD_HPP

#include <iostream>

#include "ParticleContainerBase.hpp"

#include <map>
#include <boost/lexical_cast.hpp>
#include <boost/array.hpp>
#include <boost/utility/enable_if.hpp>
#include <boost/type_traits/is_same.hpp>
#include <boost/scoped_ptr.hpp>
#include <boost/foreach.hpp>
#include "exceptions.hpp"
#include "generator.hpp"
#include "filters.hpp"
#include "Particle.hpp"
#include "ParticleID.hpp"
#include "SpeciesTypeID.hpp"
#include "SpeciesInfo.hpp"
#include "SerialIDGenerator.hpp"
#include "Transaction.hpp"
#include "Structure.hpp"
#include "StructureType.hpp"
#include "StructureID.hpp"
#include "Surface.hpp"
#include "Region.hpp"
#include "geometry.hpp"
#include "GSLRandomNumberGenerator.hpp"
#include "Point.hpp" // XXX: workaround. should be removed later.
#include "utils/pair.hpp"
#include "utils/range.hpp"

template<typename Tderived_, typename Tlen_, typename TD_>
struct WorldTraitsBase
// Just defines some basic properties of a world. Basically just types and one constant.
// Tderived_ is the class/structure that inherits from this struct (subclass) -> WorldTraits and CyclicWorldTraits
{
    typedef std::size_t                                     size_type;
    typedef Tlen_                                           length_type;
    typedef TD_                                             D_type;
    typedef TD_                                             v_type;
    typedef ParticleID                                      particle_id_type;       // identifier type for particles
    typedef SerialIDGenerator<particle_id_type>             particle_id_generator;
    typedef SpeciesTypeID                                   species_id_type;        // identifier type for species and (structure types)
    typedef SpeciesTypeID                                   structure_type_id_type;
    typedef StructureType                                   structure_type_type;    // definition of what class to use for structure_type
//    typedef Particle<length_type, D_type, species_id_type>  particle_type;      // type for particles, NOTE why is there no v_type here?
//    typedef std::string                                     structure_id_type;  // identifier type for structures
    typedef StructureID                                     structure_id_type;
    typedef SerialIDGenerator<structure_id_type>            structure_id_generator;
    typedef Particle<length_type, D_type, species_id_type,
                     structure_id_type>                     particle_type;      // type for particles, NOTE why is there no v_type here?

    typedef SpeciesInfo<species_id_type, D_type,
                        length_type, structure_type_id_type>species_type;       // This is the definition of the species in the presence of a World
    typedef Vector3<length_type>                                                    point_type;
    typedef typename particle_type::shape_type::position_type                       position_type;
    typedef GSLRandomNumberGenerator                                                rng_type;
    typedef Structure<Tderived_>                                                    structure_type; // The structure_type is parameterized with the subclass??

    static const Real TOLERANCE = 1e-7;         // tolerance of what?
};



template<typename Tlen_, typename TD_>
struct WorldTraits: public WorldTraitsBase<WorldTraits<Tlen_, TD_>, Tlen_, TD_>
{
// This defines a structure (like a class) with some rudimentairy properties of the world
// It inherits from WorldTraitsBase which is parameterized with the fully defined WorldTraits class/structure.
//
// Implements mostly boundaryconditions related stuff (distance, apply_boundary, cyclic_transpose)

public:
    // This is the normal world (without periodic/cyclic boundary conditions)
    typedef WorldTraitsBase<WorldTraits<Tlen_, TD_>, Tlen_, TD_> base_type;
    typedef typename base_type::length_type length_type;
    typedef typename base_type::position_type position_type;

    template<typename Tval_>
    static Tval_ apply_boundary(Tval_ const& v, length_type const& world_size)
    // apply_boundary is here defined for consistency
    {
        return v;
    }

    template<typename Tval_>
    static Tval_ cyclic_transpose(Tval_ const& p0, Tval_ const& p1, length_type const& world_size)
    // ditto
    {
        return p0;
    }

    template<typename T1_, typename T2_>
    static length_type distance(T1_ const& p0, T2_ const& p1, length_type const& world_size)
    // The distance function
    {
        return ::distance(p0, p1);
    }

    template<typename Toc_, typename Tfun_, typename Tsphere_>
    static void each_neighbor(Toc_& oc, Tfun_& fun, Tsphere_ const& pos)
    {
        oc.each_neighbor(oc.index(pos), fun);
    }

    template<typename Toc_, typename Tfun_, typename Tsphere_>
    static void each_neighbor(Toc_ const& oc, Tfun_& fun, Tsphere_ const& pos)
    {
        oc.each_neighbor(oc.index(pos), fun);
    }

    template<typename Toc_, typename Tfun_, typename Tsphere_>
    static void take_neighbor(Toc_& oc, Tfun_& fun, const Tsphere_& cmp)
    {
        take_neighbor(oc, fun, cmp);
    }

    template<typename Toc_, typename Tfun_, typename Tsphere_>
    static void take_neighbor(Toc_ const& oc, Tfun_& fun, const Tsphere_& cmp)
    {
        take_neighbor(oc, fun, cmp);
    }
};

template<typename Tlen_, typename TD_>
struct CyclicWorldTraits: public WorldTraitsBase<CyclicWorldTraits<Tlen_, TD_>, Tlen_, TD_>
{
// Use these world traits for a world WITH periodic/cyclic boundary conditions

public:
    // Specify the types being used.
    typedef WorldTraitsBase<CyclicWorldTraits<Tlen_, TD_>, Tlen_, TD_> base_type;
    typedef typename base_type::length_type                            length_type;
    typedef typename base_type::position_type                          position_type;

    template<typename Tval_>
    static Tval_ apply_boundary(Tval_ const& v, length_type const& world_size)
    // This applied the periodic boundary conditions, unclear what kind of type 'Tval_' is.
    {
        return ::apply_boundary(v, world_size);
    }

    static length_type cyclic_transpose(length_type const& p0, length_type const& p1, length_type const& world_size)
    // selects the copy of p0 (over the periodic boundaries) such that the distance between p0 and p1 is minimized over
    // over the periodic boundary conditions.
    // This is the length_type version
    {
        return ::cyclic_transpose(p0, p1, world_size);
    }

    static position_type cyclic_transpose(position_type const& p0, position_type const& p1, length_type const& world_size)
    // This is the position_type version
    {
        return ::cyclic_transpose(p0, p1, world_size);
    }

    template<typename T1_, typename T2_>
    static length_type distance(T1_ const& p0, T2_ const& p1, length_type const& world_size)
    // This distance function in the presence of the periodic boundary conditions
    {
        return distance_cyclic(p0, p1, world_size);
    }

    template<typename Toc_, typename Tfun_, typename Tsphere_>
    static void each_neighbor(Toc_& oc, Tfun_& fun, Tsphere_ const& pos)
    {
        oc.each_neighbor_cyclic(oc.index(pos), fun);
    }

    template<typename Toc_, typename Tfun_, typename Tsphere_>
    static void each_neighbor(Toc_ const& oc, Tfun_& fun, Tsphere_ const& pos)
    {
        oc.each_neighbor_cyclic(oc.index(pos), fun);
    }

    template<typename Toc_, typename Tfun_, typename Tsphere_>
    static void take_neighbor(Toc_& oc, Tfun_& fun, const Tsphere_& cmp)
    {
        take_neighbor_cyclic(oc, fun, cmp);
    }

    template<typename Toc_, typename Tfun_, typename Tsphere_>
    static void take_neighbor(Toc_ const& oc, Tfun_& fun, const Tsphere_& cmp)
    {
        take_neighbor_cyclic(oc, fun, cmp);
    }
};

template<typename Ttraits_>
class World: public ParticleContainerBase<World<Ttraits_>, Ttraits_>
{
// The world 'is a' ParticleContainerBase which is parameterize with the same World class
// (I think the concrete implementation of it, not just the template)
// The world is a template that is parametrized with the traits (these define a number of datatypes used in the definition of
// the world?)
public:
    typedef Ttraits_                        traits_type;                // Needs to be of WorldTraits type???
    typedef ParticleContainerBase<World>    base_type;                  // the class from which is inherited -> base class
    typedef ParticleContainer<traits_type>  particle_container_type;    // ParticleContainer is again the superclass of ParticleContainerBase

    // define some short hand names for all the properties (traits) of the world
    typedef typename traits_type::length_type               length_type;
    typedef typename traits_type::species_type              species_type;
    typedef typename traits_type::position_type             position_type;
    typedef typename traits_type::particle_type             particle_type;
    typedef typename traits_type::particle_id_type          particle_id_type;
    typedef typename traits_type::particle_id_generator     particle_id_generator;
    typedef typename traits_type::species_id_type           species_id_type;        // Note that this is also the structure_type_id_type
    typedef typename traits_type::particle_type::shape_type particle_shape_type;
    typedef typename traits_type::size_type                 size_type;
    typedef typename traits_type::structure_id_type         structure_id_type;
    typedef typename traits_type::structure_id_generator    structure_id_generator;
    typedef typename traits_type::structure_type            structure_type;
    typedef typename traits_type::structure_type_type       structure_type_type;
    typedef typename traits_type::structure_type_id_type    structure_type_id_type;

    //
    typedef std::pair<const particle_id_type, particle_type>                    particle_id_pair;      // defines the pid_particle_pair tuple
    typedef std::pair<const structure_id_type, structure_type>                  structure_id_pair;
    typedef std::pair<position_type, length_type>                               projected_type;
    typedef typename particle_container_type::structure_id_and_distance_pair    structure_id_and_distance_pair;

protected:
    typedef std::map<species_id_type, species_type>                         species_map;
    typedef std::set<particle_id_type>                                      particle_id_set;
    typedef std::map<species_id_type, particle_id_set>                      per_species_particle_id_set;
    typedef select_second<typename species_map::value_type>                 species_second_selector_type;   // not sure what this does.

    typedef std::map<structure_type_id_type, structure_type_type>           structure_type_map;
    typedef std::map<structure_id_type, boost::shared_ptr<structure_type> > structure_map;          // note: this is a structure_map_type
    typedef std::set<structure_id_type>                                     structure_id_set;
    typedef std::map<structure_type_id_type, structure_id_set>              per_structure_type_structure_id_set;
    typedef select_second<typename structure_map::value_type>               structures_second_selector_type;
    typedef select_second<typename structure_type_map::value_type>          structure_types_second_selector_type;

public:
    typedef boost::transform_iterator<species_second_selector_type,
            typename species_map::const_iterator>                           species_iterator;
    typedef boost::transform_iterator<structures_second_selector_type,
            typename structure_map::const_iterator>                         structure_iterator;
    typedef boost::transform_iterator<structure_types_second_selector_type,
            typename structure_type_map::const_iterator>                    structure_type_iterator;
    typedef sized_iterator_range<species_iterator>                          species_range;
    typedef sized_iterator_range<structure_iterator>                        structures_range;
    typedef sized_iterator_range<structure_type_iterator>                   structure_types_range;

public:
    // The constructor
    World(length_type world_size = 1., size_type size = 1)
        : base_type(world_size, size)
    {
        default_structure_id_ = structidgen_();
    }
    // TODO Add the default structure of the default structure_type here?

    // To create new particles
    virtual particle_id_pair new_particle(species_id_type const& sid,
            position_type const& pos)
    {
        species_type const& species(get_species(sid));
        particle_id_pair retval(pidgen_(),
                                particle_type(sid, particle_shape_type(pos, species.radius()),
                                              default_structure_id_, species.D(), species.v() ));

        update_particle(retval);
        return retval;
    }
    // To update particles
    virtual bool update_particle(particle_id_pair const& pi_pair)
    {
        typename base_type::particle_matrix_type::iterator i(
                base_type::pmat_.find(pi_pair.first));
        if (i != base_type::pmat_.end())
        {
            if ((*i).second.sid() != pi_pair.second.sid())
            // If the species changed we need to erase it from the 'species_id -> particle ids' mapping
            {
                particle_pool_[(*i).second.sid()].erase((*i).first);
                particle_pool_[pi_pair.second.sid()].insert(pi_pair.first);
            }
            base_type::pmat_.update(i, pi_pair);
            return false;
        }
        BOOST_ASSERT(base_type::update_particle(pi_pair));
        particle_pool_[pi_pair.second.sid()].insert(pi_pair.first);
        return true;
    }
    // Remove the particle by id 'id' if present in the world.
    virtual bool remove_particle(particle_id_type const& id)
    {
        bool found(false);
        particle_id_pair pp(get_particle(id, found));
        if (!found)
        {
            return false;
        }
        particle_pool_[pp.second.sid()].erase(id);
        base_type::remove_particle(id);
        return true;
    }

    // Adds a species to the world.
    void add_species(species_type const& species)
    {
        species_map_[species.id()] = species;
        particle_pool_[species.id()] = particle_id_set();
    }
    // Get the species in the world by id
    virtual species_type const& get_species(species_id_type const& id) const
    {
        typename species_map::const_iterator i(species_map_.find(id));
        if (species_map_.end() == i)
        {
            throw not_found(std::string("Unknown species (id=") + boost::lexical_cast<std::string>(id) + ")");
        }
        return (*i).second;
    }
    // Get all the species
    species_range get_species() const
    {
        return species_range(
            species_iterator(species_map_.begin(), species_second_selector_type()),
            species_iterator(species_map_.end(), species_second_selector_type()),
            species_map_.size());
    }

    // TODO
    // remove structure
    // Add a structure 
    bool add_structure(boost::shared_ptr<structure_type> structure)
    {
//        return structure_map_.insert(std::make_pair(structure->id(), structure)).second;
        // FIXME check that the structure_type that is defined in the structure exists!
        // TODO add mapping from structure id -> set of particles
        return update_structure(std::make_pair(structidgen_(), structure))
    }


    // update structure
    // FIXME what to do when the structure has particles and moves or changes structure_type?!?!
    virtual bool update_structure(structure_id_pair const& structid_pair)
    {
        typename structure_map::const_iterator i(structure_map_.find(structid_pair.first));
        if (i != structure_map_.end())
        // The item was already found in the map
        {
            if ((*i).second.sid() != structid_pair.second.sid())
            // If the structuretype changed we need to update the 'structure_type_id -> structure ids' mapping
            {
                structure_pool_[(*i).second.sid()].erase((*i).first);
                structure_pool_[structid_pair.second.sid()].insert(structid_pair.first);
            }
            structure_map_.update(i, structid_pair);
            return false;
        }

        // create a new item in the structure mapping
        structure_map_[structid_pair.first] = structid_pair.second;
        // update the mapping 'structure_type_id -> set of structure ids'
        structure_pool_[structid_pair.second.sid()].insert(structid_pair.first);
        return true;
    }
    // Get a structure by id
    virtual boost::shared_ptr<structure_type> get_structure(structure_id_type const& id) const
    {
        typename structure_map::const_iterator i(structure_map_.find(id));
        if (structure_map_.end() == i)
        {
            throw not_found(std::string("Unknown structure (id=") + boost::lexical_cast<std::string>(id) + ")");
        }
        return (*i).second;
    }
    // Get all the structures
    virtual structures_range get_structures() const
    {
        return structures_range(
            structure_iterator(structure_map_.begin(), structures_second_selector_type()),
            structure_iterator(structure_map_.end(),   structures_second_selector_type()),
            structure_map_.size());
    }
    // Get all the structure ids by structure_type id
    structure_id_set get_structure_ids(structure_type_id_type const& sid) const
    {
        typename per_structure_type_structure_id_set::const_iterator i(
            structure_pool_.find(sid));
        if (i == structure_pool_.end())
        {
            throw not_found(std::string("Unknown structure_type (id=") + boost::lexical_cast<std::string>(sid) + ")");
        }
        return (*i).second;
    }



    virtual structure_type_id_type get_def_structure_type_id() const
    {
        return default_structure_type_id_;
    }

    virtual structure_type_type get_structure_type(structure_type_id_type const& sid) const
    {
        return default_structure_; // TODO
    }

    virtual structure_id_type get_def_structure_id() const
    {   return default_structure_id_;
    }

    // TODO
    // Add structureType
    void add_structure_type(structure_type_type const& structure_type)
    {
        structure_type_map_[structure_type.id()] = structure_type;
        structure_pool_[structure_type.id()] = structure_id_set();
    }
    // Get structureType by id
    virtual structure_type_type const& get_structure_type(structure_type_id_type const& id) const
    {
        typename structure_type_map::const_iterator i(structure_type_map_.find(id));
        if (structure_type_map_.end() == i)
        {
            throw not_found(std::string("Unknown structure_type (id=") + boost::lexical_cast<std::string>(id) + ")");
        }
        return (*i).second;
    }
    // Get all structureTypes in the world
    structure_types_range get_structure_types() const
    {
        return structure_types_range(
            structure_type_iterator(structure_type_map_.begin(), structure_types_second_selector_type()),
            structure_type_iterator(structure_type_map_.end(), structure_types_second_selector_type()),
            structure_type_map_.size());
    }

    // Get the closest surface
    virtual structure_id_and_distance_pair get_closest_surface(position_type const& pos, structure_id_type const& ignore) const
    {       
        length_type const world_size( base_type::world_size() );
        length_type ret_dist( std::numeric_limits<length_type>::max() );
        structure_id_type ret_id( ignore );

        BOOST_FOREACH(boost::shared_ptr<structure_type> const structure, get_structures())
        {
            if( structure->id() == ignore )
                continue;
            
            position_type const cyc_pos(cyclic_transpose(pos, structure->structure_position(), world_size));
            length_type const dist( fabs( structure->projected_point_on_surface( cyc_pos ).second ) );
                
            if( dist < ret_dist )
            {
                ret_dist = dist;
                ret_id = structure->id();
            }  
        }
        
        return structure_id_and_distance_pair( ret_id , ret_dist );
    }
   
    // Get all the particle ids by species
    particle_id_set get_particle_ids(species_id_type const& sid) const
    {
        typename per_species_particle_id_set::const_iterator i(
            particle_pool_.find(sid));
        if (i == particle_pool_.end())
        {
            throw not_found(std::string("Unknown species (id=") + boost::lexical_cast<std::string>(sid) + ")");
        }
        return (*i).second;
    }

///////////// Member variables
private:
    Structure_id_type           default_structure_id_;
    structure_type_id_type      default_structure_type_id_;
    structure_type_type         default_structure_; // TODO remove later again.
    particle_id_generator               pidgen_;            // generator used to produce the unique ids for the particles
    structure_id_generator              structidgen_;
    species_map                         species_map_;       // mapping: species_id -> species
    structure_type_map                  structure_type_map_;// mapping: structure_type_id -> structure_type
    structure_map                       structure_map_;     // mapping: structure_id -> structure
    per_structure_type_structure_id_set structure_pool_;    // mapping: structure_type_id -> set of structure ids of that structure type
    per_species_particle_id_set         particle_pool_;     // mapping: species_id -> set of particle ids of that species
};

#endif /* WORLD_HPP */
