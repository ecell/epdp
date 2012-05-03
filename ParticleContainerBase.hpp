#ifndef PARTICLE_CONTAINER_BASE_HPP
#define PARTICLE_CONTAINER_BASE_HPP
#include "utils/range.hpp"
#include "utils/get_mapper_mf.hpp"
#include "utils/unassignable_adapter.hpp"
#include "MatrixSpace.hpp"
#include "abstract_set.hpp"
#include "generator.hpp"
#include "exceptions.hpp"
#include "ParticleContainer.hpp"
#include "StructureContainer.hpp"
#include "Transaction.hpp"

template<typename Ttraits_>
struct ParticleContainerUtils
{
    // The structure is parameterized with the traits of the world
    typedef Ttraits_ traits_type;
    // shorthand names for all the types that we use
    typedef typename traits_type::length_type                   length_type;
    typedef typename traits_type::particle_type                 particle_type;
    typedef typename traits_type::particle_id_type              particle_id_type;
    typedef std::pair<const particle_id_type, particle_type>    particle_id_pair;
    typedef std::pair<particle_id_pair, length_type>            particle_id_pair_and_distance;

    template<typename Tobject_and_distance_list_>
    struct distance_comparator:
            public std::binary_function<
                typename Tobject_and_distance_list_::placeholder,
                typename Tobject_and_distance_list_::placeholder,
                bool>
    {
        typedef Tobject_and_distance_list_              list_type;
        typedef typename list_type::placeholder         first_argument_type;
        typedef typename list_type::const_caster        const_caster;
        bool operator()(first_argument_type const& lhs,
                        first_argument_type const& rhs) const
        {
            return c_(lhs).second < c_(rhs).second;
        }

        const_caster c_;
    };



    template<typename Tobject_and_distance_list_, typename Tset_>
    struct overlap_checker
    {
        typedef Tobject_and_distance_list_  list_type;

        // The constructor
        overlap_checker(Tset_ const& ignore = Tset_()): ignore_(ignore), result_(0) {}

        template<typename Titer_>
        void operator()(Titer_ const& i, length_type const& dist)
        {
            if (!contains(ignore_, (*i).first))
            {
                if (!result_)
                {
                    result_ = new list_type();
                }
                result_->push_back(std::make_pair(*i, dist));
            }
        }

        // The list of items is only sorted when the results are queried.
        list_type* result() const
        {
            if (result_)
            {
                std::sort(result_->pbegin(), result_->pend(), compare_);
            }
            return result_;
        }

    private:
        Tset_ const&                    ignore_;
        list_type*                      result_;
        distance_comparator<list_type>  compare_;
    };
};


template<typename Tderived_, typename Ttraits_ = typename Tderived_::traits_type>
class ParticleContainerBase
    : public ParticleContainer<Ttraits_>
{
// This inherits from the ParticleContainer class which is just an abstract data type.
// Here most of the methods of the ParticleContainer are actually implemented.
public:
    typedef ParticleContainerUtils<Ttraits_> utils;
    typedef ParticleContainer<Ttraits_> base_type;
    typedef Ttraits_ traits_type;

    // define some shorthands for all the types from the traits that we use.
    typedef typename traits_type::length_type               length_type;
    typedef typename traits_type::species_type              species_type;
    typedef typename traits_type::position_type             position_type;
    typedef typename traits_type::particle_type             particle_type;
    typedef typename traits_type::particle_id_type          particle_id_type;
    typedef typename traits_type::particle_id_generator     particle_id_generator;
    typedef typename traits_type::species_id_type           species_id_type;
    typedef typename traits_type::particle_type::shape_type particle_shape_type;
    typedef typename traits_type::size_type                 size_type;
    typedef typename traits_type::structure_id_type         structure_id_type;
    typedef typename traits_type::structure_type            structure_type;
    typedef typename base_type::particle_id_set             particle_id_set;
    typedef typename base_type::particle_id_pair            particle_id_pair;
    typedef typename base_type::structure_id_set            structure_id_set;
    typedef typename base_type::structure_id_pair           structure_id_pair;
    typedef typename base_type::structure_types_range       structure_types_range;

    typedef typename base_type::structure_map               structure_map;
    typedef typename base_type::structures_range            structures_range;
    typedef Transaction<traits_type>                        transaction_type;

    typedef MatrixSpace<particle_type, particle_id_type, get_mapper_mf> particle_matrix_type;
    typedef abstract_limited_generator<particle_id_pair>                particle_id_pair_generator;
    typedef std::pair<particle_id_pair, length_type>                    particle_id_pair_and_distance;
    typedef sized_iterator_range<typename particle_matrix_type::const_iterator> particle_id_pair_range;

    typedef StructureContainer<structure_type, structure_id_type>       structure_container_type;

    typedef typename base_type::particle_id_pair_and_distance_list                              particle_id_pair_and_distance_list;
    typedef typename base_type::structure_id_pair_and_distance_list                             structure_id_pair_and_distance_list;
    typedef typename base_type::structure_id_pair_and_distance                                  structure_id_pair_and_distance;


public:
//    ParticleContainerBase(length_type world_size, size_type size)
//        : pmat_(world_size, size) {}
    ParticleContainerBase(length_type world_size, size_type size, structure_id_type default_structid)
        : pmat_(world_size, size), structures_(default_structid) {}
    // constructors

    virtual size_type num_particles() const
    {
        return pmat_.size();
    }

    virtual length_type world_size() const
    {
        return pmat_.world_size();
    }

    length_type cell_size() const
    {
        return pmat_.cell_size();
    }

    size_type matrix_size() const
    {
        return pmat_.matrix_size();
    }

    template<typename T_>
    length_type distance(T_ const& lhs, position_type const& rhs) const
    {
        return traits_type::distance(lhs, rhs, world_size());
    }

    virtual length_type distance(position_type const& lhs,
                                 position_type const& rhs) const
    {
        return traits_type::distance(lhs, rhs, world_size());
    }

    virtual position_type apply_boundary(position_type const& v) const
    {
        return traits_type::apply_boundary(v, world_size());
    }

    virtual length_type apply_boundary(length_type const& v) const
    {
        return traits_type::apply_boundary(v, world_size());
    }

    virtual position_type cyclic_transpose(position_type const& p0, position_type const& p1) const
    {
        return traits_type::cyclic_transpose(p0, p1, world_size());
    }

    virtual length_type cyclic_transpose(length_type const& p0, length_type const& p1) const
    {
        return traits_type::cyclic_transpose(p0, p1, world_size());
    }

    // THIS SEEMS STRANGE TO PUT THIS HERE. 
    template<typename T1_>
    T1_ calculate_pair_CoM(
        T1_ const& p1, T1_ const& p2, 
        typename element_type_of<T1_>::type const& D1,
        typename element_type_of<T1_>::type const& D2)
    {
        typedef typename element_type_of< T1_ >::type element_type;   

        T1_ retval;

        const T1_ p2t(cyclic_transpose(p2, p1));

        return modulo(
            divide(
                add(multiply(p1, D2), multiply(p2t, D1)),
                add(D1, D2)),
            world_size());
    }

    virtual particle_id_pair_and_distance_list* check_overlap(particle_shape_type const& s) const
    {
        return check_overlap<particle_shape_type>(s);
    }

    virtual particle_id_pair_and_distance_list* check_overlap(particle_shape_type const& s, particle_id_type const& ignore) const
    {
        return check_overlap(s, array_gen(ignore)); // Why no parameterization of the template here?
    }

    virtual particle_id_pair_and_distance_list* check_overlap(particle_shape_type const& s, particle_id_type const& ignore1, particle_id_type const& ignore2) const
    {
        return check_overlap(s, array_gen(ignore1, ignore2));
    }

    template<typename Tsph_, typename Tset_>
    particle_id_pair_and_distance_list* check_overlap(Tsph_ const& s, Tset_ const& ignore,
        typename boost::disable_if<boost::is_same<Tsph_, particle_id_pair> >::type* =0) const
    {
        typename utils::template overlap_checker<particle_id_pair_and_distance_list, Tset_> oc(ignore);
        traits_type::take_neighbor(pmat_, oc, s);
        return oc.result();
    }
    // This is the same method as the one above but now without the 'ignore' list
    template<typename Tsph_>
    particle_id_pair_and_distance_list* check_overlap(Tsph_ const& s,
        typename boost::disable_if<boost::is_same<Tsph_, particle_id_pair> >::type* =0) const
    {
        typename utils::template overlap_checker<particle_id_pair_and_distance_list, boost::array<particle_id_type, 0> > oc;
        traits_type::take_neighbor(pmat_, oc, s);
        return oc.result();
    }

    virtual structure_id_pair_and_distance_list* check_surface_overlap(particle_shape_type const& s, position_type const& old_pos, structure_id_type const& current,
                                                                       length_type const& sigma) const
    {
        typename utils::template overlap_checker<structure_id_pair_and_distance_list, boost::array<structure_id_type, 0> > checker;
        const structure_id_set visible_structures (structures_.get_visible_structures(current));

        structure_map temp_map;
        for (typename structure_id_set::const_iterator i(visible_structures.begin()),
                                                       e(visible_structures.end());
             i != e; ++i)
        {
            temp_map[(*i)] = get_structure(*i);
        }

        for (typename structure_map::const_iterator i(temp_map.begin()),
                                                    e(temp_map.end());
             i != e; ++i)
        {
            const position_type cyc_pos(cyclic_transpose(s.position(), ((*i).second)->position())); //
            const length_type dist((*i).second->newBD_distance(cyc_pos, s.radius(), old_pos, sigma));
            if (dist < s.radius())
            {
                checker(i, dist);
            }
        }
        return checker.result();
    }

    particle_id_pair get_particle(particle_id_type const& id, bool& found) const
    {
        typename particle_matrix_type::const_iterator i(pmat_.find(id));
        if (pmat_.end() == i) {
            found = false;
            return particle_id_pair();
        }
        found = true;
        return *i;
    }

    virtual particle_id_pair get_particle(particle_id_type const& id) const
    {
        typename particle_matrix_type::const_iterator i(pmat_.find(id));
        if (pmat_.end() == i) {
            throw not_found(std::string("No such particle: id=")
                    + boost::lexical_cast<std::string>(id));
        }
        return *i;
    }

    virtual bool has_particle(particle_id_type const& id) const
    {
        return pmat_.end() != pmat_.find(id);
    }

    virtual transaction_type* create_transaction();     // The implementation is below as an inline function?

    virtual particle_id_pair_generator* get_particles() const
    {
        return make_range_generator<particle_id_pair>(pmat_);
    }

    particle_id_pair_range get_particles_range() const
    {
        return particle_id_pair_range(pmat_.begin(), pmat_.end(), pmat_.size());
    }

    virtual bool update_particle(particle_id_pair const& pi_pair)
    {
        return pmat_.update(pi_pair).second;
    }

    virtual bool remove_particle(particle_id_type const& id)
    {
        return pmat_.erase(id);
    }

    ///// Structure stuff
    virtual bool has_structure(structure_id_type const& id) const
    {
        return structures_.has_structure(id);
    }
    virtual boost::shared_ptr<structure_type> get_structure(structure_id_type const& id) const
    {
        return structures_.get_structure(id);
    }
    virtual structures_range get_structures() const
    {
        return structures_.get_structures_range();
    }
    virtual bool update_structure(structure_id_pair const& structid_pair)
    {
        return structures_.update_structure(structid_pair);
    }
    virtual bool remove_structure(structure_id_type const& id)
    {
        return structures_.remove_structure(id);
    }
    virtual structure_id_pair_and_distance_list* get_close_structures(position_type const& pos, structure_id_type const& current_struct_id,
                                                                      structure_id_type const& ignore) const
    {
        typename utils::template overlap_checker<structure_id_pair_and_distance_list, boost::array<structure_id_type, 1> > checker(array_gen(ignore));
        const structure_id_set visible_structures (structures_.get_visible_structures(current_struct_id));

        structure_map temp_map;
        for (typename structure_id_set::const_iterator i(visible_structures.begin()),
                                                       e(visible_structures.end());
             i != e; ++i)
        {
            temp_map[(*i)] = get_structure(*i);
        }

        for (typename structure_map::const_iterator i(temp_map.begin()),
                                                    e(temp_map.end());
             i != e; ++i)
        {
            const position_type cyc_pos(cyclic_transpose(pos, ((*i).second)->position()));
            const length_type dist((*i).second->distance(cyc_pos));
            checker(i, dist);
        }
        return checker.result();
    }

///////// Member variables
protected:
    particle_matrix_type        pmat_;          // the structure (MatrixSpace) containing the particles.
    structure_container_type    structures_;    // an object containing the structures.
};


//////// Inline methods are defined separately
template<typename Tderived_, typename Ttraits_>
inline Transaction<Ttraits_>*
ParticleContainerBase<Tderived_, Ttraits_>::create_transaction()
{
    return new TransactionImpl<ParticleContainerBase>(*this);
}

#endif /* PARTICLE_CONTAINER_BASE_HPP */
