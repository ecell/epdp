#ifndef WORLD_HPP
#define WORLD_HPP

#include "ParticleContainerBase.hpp"

#include <map>
#include <boost/lexical_cast.hpp>
#include <boost/array.hpp>
#include <boost/foreach.hpp>
#include <boost/utility/enable_if.hpp>
#include <boost/type_traits/is_same.hpp>
#include <boost/type_traits/is_pointer.hpp>
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
#include "StructureID.hpp"
#include "geometry.hpp"
#include "GSLRandomNumberGenerator.hpp"
#include "Point.hpp" // XXX: workaround. should be removed later.
#include "utils/pair.hpp"

template<typename Tderived_, typename Tlen_, typename TD_>
struct WorldTraitsBase
{
    typedef std::size_t size_type;
    typedef Tlen_ length_type;
    typedef TD_ D_type;
    typedef TD_ v_type;
    typedef ParticleID particle_id_type;
    typedef SerialIDGenerator<particle_id_type> particle_id_generator_type;
    typedef SpeciesTypeID species_id_type;
    typedef Structure<Tderived_> structure_type;
    typedef StructureID structure_id_type;
    typedef SerialIDGenerator<structure_id_type> structure_id_generator_type;
    typedef std::string structure_type_id_type;
    typedef Particle<length_type, D_type, species_id_type, structure_id_type> particle_type;
    typedef SpeciesInfo<species_id_type, D_type, length_type, structure_type_id_type> species_type;
    typedef Vector3<length_type> point_type;
    typedef typename particle_type::shape_type::position_type position_type;
    typedef GSLRandomNumberGenerator rng_type;

    static const Real TOLERANCE = 1e-7;
};

template<typename Tlen_, typename TD_>
struct WorldTraits: public WorldTraitsBase<WorldTraits<Tlen_, TD_>, Tlen_, TD_>
{
public:
    typedef WorldTraitsBase<WorldTraits<Tlen_, TD_>, Tlen_, TD_> base_type;
    typedef typename base_type::length_type length_type;
    typedef typename base_type::position_type position_type;

    template<typename Tval_>
    static Tval_ apply_boundary(Tval_ const& v, length_type const& world_size)
    {
        return v;
    }

    template<typename Tval_>
    static Tval_ cyclic_transpose(Tval_ const& p0, Tval_ const& p1, length_type const& world_size)
    {
        return p0;
    }

    template<typename T1_, typename T2_>
    static length_type distance(T1_ const& p0, T2_ const& p1, length_type const& world_size)
    {
        return ::distance(p0, p1);
    }

    template<typename T1_, typename T2_>
    static length_type distance(boost::shared_ptr<T1_> const& p0, T2_ const& p1, length_type const& world_size)
    {
        return p0->distance(p1);
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
public:
    typedef WorldTraitsBase<CyclicWorldTraits<Tlen_, TD_>, Tlen_, TD_> base_type;
    typedef typename base_type::length_type length_type;
    typedef typename base_type::position_type position_type;

    template<typename Tval_>
    static Tval_ apply_boundary(Tval_ const& v, length_type const& world_size)
    {
        return ::apply_boundary(v, world_size);
    }

    static length_type cyclic_transpose(length_type const& p0, length_type const& p1, length_type const& world_size)
    {
        return ::cyclic_transpose(p0, p1, world_size);
    }

    static position_type cyclic_transpose(position_type const& p0, position_type const& p1, length_type const& world_size)
    {
        return ::cyclic_transpose(p0, p1, world_size);
    }

    template<typename T1_, typename T2_>
    static length_type distance(T1_ const& p0, T2_ const& p1, length_type const& world_size)
    {
        return distance_cyclic(p0, p1, world_size);
    }

    template<typename T1_, typename T2_>
    static length_type distance(boost::shared_ptr<T1_> const& p0, T2_ const& p1, length_type const& world_size)
    {
        return p0->distance_cyclic(p1, world_size);
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
public:
    typedef Ttraits_ traits_type;
    typedef ParticleContainerBase<World> base_type;
    typedef ParticleContainer<traits_type> particle_container_type;
    typedef typename traits_type::length_type length_type;
    typedef typename traits_type::species_type species_type;
    typedef typename traits_type::position_type position_type;
    typedef typename traits_type::particle_type particle_type;
    typedef typename traits_type::particle_id_type particle_id_type;
    typedef typename traits_type::particle_id_generator_type particle_id_generator_type;
    typedef typename traits_type::species_id_type species_id_type;
    typedef typename traits_type::particle_type::shape_type particle_shape_type;
    typedef typename traits_type::size_type size_type;
    typedef typename traits_type::structure_type_id_type structure_type_id_type;
    typedef typename traits_type::structure_type structure_type;
    typedef typename traits_type::structure_id_type structure_id_type;
    typedef typename traits_type::structure_id_generator_type structure_id_generator_type;
    typedef std::pair<const structure_id_type, boost::shared_ptr<structure_type> > structure_id_pair;
    typedef std::pair<const particle_id_type, particle_type> particle_id_pair;

    typedef abstract_limited_generator<structure_id_pair> structure_id_pair_generator;
    typedef std::pair<structure_id_pair, length_type> structure_id_pair_and_distance;
    typedef unassignable_adapter<structure_id_pair_and_distance, get_default_impl::std::vector> structure_id_pair_and_distance_list;

protected:
    typedef std::map<species_id_type, species_type> species_map;
    typedef std::set<particle_id_type> particle_id_set;
    typedef std::map<species_id_type, particle_id_set> per_species_particle_id_set;
    typedef select_second<typename species_map::value_type> species_second_selector_type;

    // typedef std::map<structure_type_id_type, boost::shared_ptr<structure_type> > structure_map;
    typedef std::map<structure_id_type, boost::shared_ptr<structure_type> > structure_map;
    typedef std::map<structure_type_id_type, structure_id_type> structure_id_map;
    typedef select_first<typename structure_map::value_type> structure_first_selector_type;
    typedef std::set<structure_id_type> structure_id_set;
    // typedef select_second<typename structure_map::value_type> structure_second_selector_type;
 
public:
    typedef boost::transform_iterator<species_second_selector_type,
            typename species_map::const_iterator> species_iterator;
    typedef sized_iterator_range<species_iterator> species_range;

    // typedef boost::transform_iterator<structure_second_selector_type,
    //         typename structure_map::const_iterator> structure_iterator;
    // typedef sized_iterator_range<structure_iterator> structures_range;
    typedef sized_iterator_range<typename structure_map::const_iterator> structures_range;

public:
    World(length_type world_size = 1., size_type size = 1)
        : base_type(world_size, size) {}

    virtual particle_id_pair new_particle(species_id_type const& sid,
            position_type const& pos)
    {
        species_type const& species(get_species(sid));
        structure_id_pair stp(get_structure(species.structure_id()));
        particle_id_pair retval(particle_id_generator_(),
            particle_type(sid, particle_shape_type(pos, species.radius()),
                          species.D(), stp.first, species.v()));
        update_particle(retval);
        return retval;
    }

    virtual bool update_particle(particle_id_pair const& pi_pair)
    {
        typename base_type::particle_matrix_type::iterator i(
                base_type::pmat_.find(pi_pair.first));
        if (i != base_type::pmat_.end())
        {
            if ((*i).second.sid() != pi_pair.second.sid())
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

    void add_species(species_type const& species)
    {
        species_map_[species.id()] = species;
        particle_pool_[species.id()] = particle_id_set();
    }

    virtual species_type const& get_species(species_id_type const& id) const
    {
        typename species_map::const_iterator i(species_map_.find(id));
        if (species_map_.end() == i)
        {
            throw not_found(std::string("Unknown species (id=") + boost::lexical_cast<std::string>(id) + ")");
        }
        return (*i).second;
    }

    species_range get_species() const
    {
        return species_range(
            species_iterator(species_map_.begin(), species_second_selector_type()),
            species_iterator(species_map_.end(), species_second_selector_type()),
            species_map_.size());
    }

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

    structure_id_pair add_structure(boost::shared_ptr<structure_type> structure)
    {
        structure_id_type structure_id(structure_id_generator_());
        structure_id_map_.insert(std::make_pair(structure->id(), structure_id));
        return (*structure_map_.insert(
                    std::make_pair(structure_id, structure)).first);
    }

    virtual structure_id_pair get_structure(structure_id_type const& id) const
    {
        typename structure_map::const_iterator i(structure_map_.find(id));
        if (structure_map_.end() == i)
        {
            throw not_found(std::string("Unknown structure (id=") + boost::lexical_cast<std::string>(id) + ")");
        }
        return (*i);
    }

    virtual structure_id_pair get_structure(structure_type_id_type const& id) const
    {
        typename structure_id_map::const_iterator i(
            structure_id_map_.find(id));
        if (structure_id_map_.end() == i)
        {
            throw not_found(std::string("Unknown structure (id=") + boost::lexical_cast<std::string>(id) + ")");
        }
        return get_structure((*i).second);
    }

    structures_range get_structures() const
    {
        return structures_range(
            structure_map_.begin(), structure_map_.end(), structure_map_.size());
    }

    structure_id_set get_structure_ids() const
    {
        return structure_id_set(
            boost::make_transform_iterator(structure_map_.begin(), structure_first_selector_type()),
            boost::make_transform_iterator(structure_map_.end(), structure_first_selector_type()));
    }

    std::pair<structure_id_type, length_type> get_closest_structure(position_type const& p) const
    {
        std::pair<structure_id_type, length_type> closest(
            structure_id_type(), std::numeric_limits<length_type>::infinity());
        BOOST_FOREACH(structure_id_pair const st_pair, structure_map_)
        {
            length_type const distance_to_structure(
                fabs(base_type::distance(st_pair.second, p)));
            if (distance_to_structure < closest.second)
            {
                closest.first = st_pair.first;
                closest.second = distance_to_structure;
            }
        }
        return closest;
    }

    structure_id_pair_and_distance_list* check_overlap_with_structures(particle_shape_type const& s, structure_id_type const& ignore) const
    {
        return check_overlap_with_structures(s, array_gen(ignore));
    }

    template<typename Tsph_, typename Tset_>
    structure_id_pair_and_distance_list* check_overlap_with_structures(Tsph_ const& s, Tset_ const& ignore, typename boost::disable_if<boost::is_same<Tsph_, particle_id_pair> >::type* = 0) const
    {
        structure_id_pair_and_distance_list* result;
        result = new structure_id_pair_and_distance_list();
        BOOST_FOREACH(structure_id_pair const st_pair, structure_map_)
        {
            if (!contains(ignore, st_pair.first))
            {
                length_type const distance_to_structure(
                    fabs(base_type::distance(st_pair.second, s.position())));
                if (distance_to_structure < s.radius())
                {
                    result->push_back(
                        std::make_pair(st_pair, distance_to_structure));
                }
            }
        }
        return result;
    }
    

private:
    particle_id_generator_type particle_id_generator_;
    structure_id_generator_type structure_id_generator_;
    species_map species_map_;
    structure_map structure_map_;
    structure_id_map structure_id_map_;
    per_species_particle_id_set particle_pool_;
};

#endif /* WORLD_HPP */
