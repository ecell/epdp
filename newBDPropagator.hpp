#ifndef NEW_BD_PROPAGATOR_HPP
#define NEW_BD_PROPAGATOR_HPP

#include <algorithm>
#include <boost/bind.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/range/size.hpp>
#include <boost/range/begin.hpp>
#include <boost/range/end.hpp>
#include <boost/range/const_iterator.hpp>
#include <boost/utility/enable_if.hpp>
#include <boost/scoped_ptr.hpp>
#include <boost/shared_ptr.hpp>
#include "Defs.hpp"
#include "generator.hpp"
#include "exceptions.hpp"
#include "freeFunctions.hpp"
#include "utils.hpp"
#include "utils/random.hpp"
#include "utils/get_default_impl.hpp"
#include "Logger.hpp"

#include <iostream>

template<typename Ttraits_>
class newBDPropagator
{
public:
    typedef Ttraits_ traits_type;
    typedef typename Ttraits_::world_type::particle_container_type particle_container_type;
    typedef typename particle_container_type::species_id_type species_id_type;
    typedef typename particle_container_type::position_type position_type;
    typedef typename particle_container_type::particle_shape_type particle_shape_type;
    typedef typename particle_container_type::species_type species_type;
    typedef typename particle_container_type::length_type length_type;
    typedef typename particle_container_type::particle_id_type particle_id_type;
    typedef typename particle_container_type::particle_type particle_type;
    typedef typename particle_container_type::particle_id_pair particle_id_pair;
    typedef typename particle_container_type::structure_id_and_distance_pair structure_id_and_distance_pair;
    typedef std::vector<particle_id_type> particle_id_vector_type;
    typedef typename particle_container_type::particle_id_pair_generator particle_id_pair_generator;
    typedef typename particle_container_type::particle_id_pair_and_distance particle_id_pair_and_distance;
    typedef typename particle_container_type::particle_id_pair_and_distance_list particle_id_pair_and_distance_list;
    typedef typename particle_container_type::structure_type structure_type;
    typedef typename particle_container_type::structure_id_type structure_id_type;    
    typedef typename traits_type::world_type::traits_type::rng_type rng_type;
    typedef typename traits_type::time_type time_type;
    typedef typename traits_type::network_rules_type network_rules_type;
    typedef typename network_rules_type::reaction_rules reaction_rules;
    typedef typename network_rules_type::reaction_rule_type reaction_rule_type;
    typedef typename traits_type::reaction_record_type reaction_record_type;
    typedef typename traits_type::reaction_recorder_type reaction_recorder_type;
    typedef typename traits_type::volume_clearer_type volume_clearer_type;


public:
    template<typename Trange_>
    newBDPropagator(
        particle_container_type& tx, network_rules_type const& rules,
        rng_type& rng, time_type dt, int max_retry_count, length_type reaction_length,
        reaction_recorder_type* rrec, volume_clearer_type* vc,
        Trange_ const& particle_ids)
        : tx_(tx), rules_(rules), rng_(rng), dt_(dt),
          max_retry_count_(max_retry_count), rrec_(rrec), vc_(vc),
          queue_(), rejected_move_count_(0), reaction_length_( reaction_length )
    {
        call_with_size_if_randomly_accessible(
            boost::bind(&particle_id_vector_type::reserve, &queue_, _1),
            particle_ids);
        for (typename boost::range_const_iterator<Trange_>::type
                i(boost::begin(particle_ids)),
                e(boost::end(particle_ids)); i != e; ++i)
        {
            queue_.push_back(*i);
        }
        shuffle(rng, queue_);
    }

    bool operator()()
    {
        if (queue_.empty())
            return false;

        particle_id_type pid(queue_.back());
        queue_.pop_back();
        particle_id_pair pp(tx_.get_particle(pid));

        LOG_DEBUG(("propagating particle %s", boost::lexical_cast<std::string>(pp.first).c_str()));

        /* Consider decay reactions first. Particles can move after decay, but this is done
           inside the 'attempt_reaction' function; for this particle can only react once . */
        try
        {
            if (attempt_reaction(pp))
                return true;
        }
        catch (propagation_error const& reason)
        {
            log_.info("first-order reaction rejected (reason: %s)", reason.what());
            ++rejected_move_count_;
            return true;
        }
        
        species_type const species(tx_.get_species(pp.second.sid()));
        length_type const r0( species.radius() );
        boost::shared_ptr<structure_type> const pp_structure( tx_.get_structure( species.structure_id() ) );
        
        structure_id_and_distance_pair struct_id_distance_pair;
        length_type core_surface_distance = 0;
        position_type new_pos;
        
       /* Sample a potential move, and check if the particle _core_ has overlapped with another 
          particle or surface. If this is the case the particlem bounces, and is returned
          to it's original position. For a particle with D = 0, new_pos = old_pos */        
        if (species.D() != 0.)
        {

            // Make a step, dependent on the surface the particle lives on.
            position_type const displacement( tx_.get_structure(species.structure_id())->
                                            bd_displacement(species.v() * dt_, std::sqrt(2.0 * species.D() * dt_), rng_) );

            new_pos = tx_.apply_boundary( add( pp.second.position(), displacement ) );
            
        }
        else
        {
            new_pos = pp.second.position();
        }
        
        /* Use a spherical shape with radius = particle_radius + reaction_length.
	       We use it to check for overlaps with other particels. */
        boost::scoped_ptr<particle_id_pair_and_distance_list> overlapped(
                tx_.check_overlap(particle_shape_type( new_pos, r0 + reaction_length_ ), pp.first));
        int particles_in_overlap(overlapped ? overlapped->size(): 0);

        /* Check if the particle at new_pos overlaps with any particle cores. */        
        bool bounced = false;       
        for( int i = 0; i < particles_in_overlap; i++ )
        {
            if( overlapped->at(i).second < r0 )
            {
                bounced = true;
                break;
            }
        }
        
        /* If the particle has not bounced with a particle and lives in the bulk: 
           check for core overlap with a surface. */
        if(!bounced && pp_structure->id() == "world")
        {
            /* Check for overlap with nearest surface. */
            struct_id_distance_pair = tx_.get_closest_surface( new_pos );
            core_surface_distance = struct_id_distance_pair.second - r0;
        
            if(core_surface_distance < 0)
                bounced = true;
        }         
         
        /* If particle is bounced, check reaction_volume for reaction parners at old_pos. */     
        if(bounced)
        {
            new_pos = pp.second.position();
            
            boost::scoped_ptr<particle_id_pair_and_distance_list> overlapped_after_bounce( 
                    tx_.check_overlap( particle_shape_type( new_pos, r0 + reaction_length_ ), 
                            pp.first) );
                            
            overlapped.swap( overlapped_after_bounce );
                                            
            particles_in_overlap = overlapped ? overlapped->size(): 0; 
                
            struct_id_distance_pair = tx_.get_closest_surface( new_pos );
            core_surface_distance = struct_id_distance_pair.second - r0;  
        }
        
        /* If their is only a surface in the reaction volume, attempt interaction. */
        if(particles_in_overlap == 0 && core_surface_distance < reaction_length_ && pp_structure->id() == "world")
        {
            try
            {
                if( false /*attempt_reaction(pp, pp_surface)*/ )
                    return true;
                else                
                {
                    LOG_DEBUG(("Particle attempted an interaction with the nonreactive surface %s.", 
                        boost::lexical_cast<std::string>(pp_structure->id()).c_str()));
                    ++rejected_move_count_;
                }
            }
            catch (propagation_error const& reason)
            {
                log_.info("surface interaction rejected (reason: %s)", reason.what());
                ++rejected_move_count_;
            }
        }
         
        switch (particles_in_overlap)
        {
        case 0:
            break;

        case 1:
        {
            
        
            particle_id_pair_and_distance const& closest(overlapped->at(0));
 		        
    		try
            {
                if( attempt_reaction(pp, closest.first) )
                    return true;
                else                
                {
                    LOG_DEBUG(("Particle attempted a reaction with a nonreactive particle %s.", 
                        boost::lexical_cast<std::string>(closest.first.first).c_str()));
                    ++rejected_move_count_;
                }
            }
            catch (propagation_error const& reason)
            {
                log_.info("second-order reaction rejected (reason: %s)", reason.what());
                ++rejected_move_count_;
            }
            /* If no reaction has occurred */    
            break;
        }
        default:
            log_.info("Multiple objects (surface and/or particle) in reaction volume; moved but no reaction initiated.");
            ++rejected_move_count_;
            break;
        }
                       
        /* If the particle did not bounce, update to new position. */
        if(!bounced)
        {
                
            particle_id_pair particle_to_update( pp.first, 
                        particle_type(species.id(), particle_shape_type(new_pos, r0), species.D()) );
                
            if (vc_)
            {
                if (!(*vc_)(particle_to_update.second.shape(), particle_to_update.first))
                {
                    log_.info("propagation move rejected.");
                    return true;
                }
            }
                
            tx_.update_particle(particle_to_update);
        }

        return true;
        
    }

    std::size_t get_rejected_move_count() const
    {
        return rejected_move_count_;
    }

private:

    position_type const make_move(species_type const& species, position_type const& np, particle_id_type const& ignore)
    {
        boost::shared_ptr<structure_type> pp_structure( tx_.get_structure( species.structure_id() ) );
    
        position_type displacement( pp_structure->
                                        bd_displacement(species.v() * dt_, std::sqrt(2.0 * species.D() * dt_), rng_) );

        position_type new_pos(tx_.apply_boundary( add( np, displacement ) ));
                        
        boost::scoped_ptr<particle_id_pair_and_distance_list> overlapped(
                                    tx_.check_overlap(
                                    particle_shape_type(new_pos, species.radius()),
                                    ignore));                            
                                    
        if ( (overlapped && overlapped->size() > 0) )
            return np;
        
        structure_id_and_distance_pair const struct_id_distance_pair = tx_.get_closest_surface( new_pos );
        length_type const core_surface_distance = struct_id_distance_pair.second - s0.radius();
                                    
        if( core_surface_distance < 0 )
            return np;
            
        return new_pos;
    }

    bool attempt_reaction(particle_id_pair const& pp)
    {
        reaction_rules const& rules(rules_.query_reaction_rule(pp.second.sid()));
        if (::size(rules) == 0)
        {
            return false;
        }

        const Real rnd(rng_() / dt_);
        Real prob = 0.;

        for (typename boost::range_const_iterator<reaction_rules>::type
                i(boost::begin(rules)), e(boost::end(rules)); i != e; ++i)
        {
            reaction_rule_type const& r(*i);
            prob += r.k();
            if (prob > rnd)
            {
                typename reaction_rule_type::species_id_range products(
                        r.get_products());
                switch (::size(products))
                {
                case 0:
                    remove_particle(pp.first);
                    break;

                case 1:
                    {
                        const species_type s0(tx_.get_species(products[0]));
                        boost::shared_ptr<structure_type> s0_struct( tx_.get_structure( s0.structure_id() ) );
                        boost::shared_ptr<structure_type> pp_struct( tx_.get_structure( tx_.get_species(pp.second.sid()).structure_id() ) );
                        
                        position_type new_pos = pp.second.position();
                        /* If the product particle does NOT live of the same structure => surface -> bulk dissociation. */
                        if( s0_struct->id() != pp_struct->id() )
                        {
                            position_type const displacement( pp_struct->surface_dissociation_vector(rng_, s0.radius(), reaction_length_ ) );

                            new_pos = tx_.apply_boundary( add( pp.second.position(), displacement ) );
                        }
                      
                        const particle_id_pair new_p(
                            pp.first, particle_type(products[0],
                                particle_shape_type(new_pos,
                                                    s0.radius()),
                                                    s0.D()));
                                                    
                        boost::scoped_ptr<particle_id_pair_and_distance_list> overlapped(tx_.check_overlap(new_p.second.shape(), new_p.first));
                        if (overlapped && overlapped->size() > 0)
                        {
                            throw propagation_error("no space");
                        }

                        if (vc_)
                        {
                            if (!(*vc_)(new_p.second.shape(), pp.first))
                            {
                                throw propagation_error("no space");
                            }
                        }

                        tx_.update_particle(new_p);

                        if (rrec_)
                        {
                            (*rrec_)(
                                reaction_record_type(
                                    r.id(), array_gen(new_p.first), pp.first));
                        }
                    }
                    break;

                case 2:
                    {
                        const species_type s0(tx_.get_species(products[0])),
                                s1(tx_.get_species(products[1]));
                        const Real D01(s0.D() + s1.D());
                        const length_type r01(s0.radius() + s1.radius());
                        int i = max_retry_count_;
                        position_type np0, np1;

                        boost::shared_ptr<structure_type> pp_struct( tx_.get_structure( tx_.get_species( pp.second.sid() ).structure_id() ) ); 
                        boost::shared_ptr<structure_type> s0_struct( tx_.get_structure( s0.structure_id() ) );
                        boost::shared_ptr<structure_type> s1_struct( tx_.get_structure( s1.structure_id() ) );
                        
                        /* If the product particles do NOT live on the same structure => surface -> surface + bulk dissociation. */
                        if( s0_struct->id() != s1_struct->id() )
                        {
                            D0_Dr = s0.D() / ( s0.D() + s1.D() );
                            D1_Dr = s1.D() / ( s0.D() + s1.D() );                        
                        
                            if(s0_struct->id() == "world")
                            {
                                position_type m( pp_struct->
                                    special_geminate_dissociation_vector( rng_, s0.radius(), s1.radius(), reaction_length_ ) );
                                
                                
                                
                                np0 = tx_.apply_boundary(pp.second.position()
                                    - m * (s0.D() / D01));
                                np1 = tx_.apply_boundary(pp.second.position()
                                    + m * (s1.D() / D01));
                            }
                            else
                            {
                                position_type m( pp_struct->
                                    special_geminate_dissociation_vector( rng_, s1.radius(), s0.radius(), reaction_length_ ) );
                                
                                np0 = tx_.apply_boundary(pp.second.position()
                                    - m * (s0.D() / D01));
                                np1 = tx_.apply_boundary(pp.second.position()
                                    + m * (s1.D() / D01));
                            }
                            


                            new_pos = tx_.apply_boundary( add( pp.second.position(), displacement ) );
                        }
                        
                        //Don't forget make_move after dissociation!
                        
                        /* Create an ipv between the products which lies within the reation volume; check for free space. */
                        for (;;)
                        {
                            if (--i < 0)
                            {
                                throw propagation_error("no space");
                            } 
               
                            position_type m( pp_struct->geminate_dissociation_vector( rng_, r01, reaction_length_ ) );
                            
                            np0 = tx_.apply_boundary(pp.second.position()
                                    - m * (s0.D() / D01));
                            np1 = tx_.apply_boundary(pp.second.position()
                                    + m * (s1.D() / D01));
                           
                            boost::scoped_ptr<particle_id_pair_and_distance_list> overlapped_s0(
                                tx_.check_overlap(
                                    particle_shape_type(np0, s0.radius()),
                                    pp.first));
                            boost::scoped_ptr<particle_id_pair_and_distance_list> overlapped_s1(
                                tx_.check_overlap(
                                    particle_shape_type(np1, s1.radius()),
                                    pp.first));
                            if (!(overlapped_s0 && overlapped_s0->size() > 0) && !(overlapped_s1 && overlapped_s1->size() > 0))
                                break;
                        }

                        /* After dissociation both particles are allowed to move. */
                        np0 = make_move(s0, np0, pp.first);
                        
                        np1 = make_move(s1, np1, pp.first);                      

			            //Check wether new positions dont overlap.

                        if (vc_)
                        {
                            if (!(*vc_)(particle_shape_type(np0, s0.radius()), pp.first) || !(*vc_)(particle_shape_type(np1, s1.radius()), pp.first))
                            {
                                throw propagation_error("no space");
                            }
                        }

                        tx_.remove_particle(pp.first);
                        const particle_id_pair
                            npp0(tx_.new_particle(s0.id(), np0)),
                            npp1(tx_.new_particle(s1.id(), np1));

                        if (rrec_)
                        {
                            (*rrec_)(
                                reaction_record_type(
                                    r.id(),
                                    array_gen(npp0.first, npp1.first),
                                    pp.first));
                        }
                    }
                    break;
                default:
                    throw not_implemented("monomolecular reactions that produce more than two products are not supported");
                }
                return true;
            }
        }
        return false;
    }

    bool attempt_reaction(particle_id_pair const& pp0, particle_id_pair const& pp1)
    {
        reaction_rules const& rules(rules_.query_reaction_rule(pp0.second.sid(), pp1.second.sid()));
        if (::size(rules) == 0)
        {
            return false;
        }

        const species_type s0(tx_.get_species(pp0.second.sid())),
                s1(tx_.get_species(pp1.second.sid()));

        const Real rnd(rng_());
        Real prob = 0;

        for (typename boost::range_const_iterator<reaction_rules>::type
                i(boost::begin(rules)), e(boost::end(rules)); i != e; ++i)
        {
            reaction_rule_type const& r(*i);

            boost::shared_ptr<structure_type> pp_structure( 
               tx_.get_structure( tx_.get_species( pp0.second.sid() ).structure_id()) );  

            Real p( r.k() * dt_ / ( pp_structure->reaction_volume( s0.radius(), s1.radius(), reaction_length_ ) ) );
	    //TODO: Check if it works.
	    p /= (2 - p);
            BOOST_ASSERT(p >= 0.);
            prob += p;
            if (prob >= 1.)
            {
                throw propagation_error(
                    "invalid acceptance ratio ("
                    + boost::lexical_cast<std::string>(p)
                    + ") for reaction rate "
                    + boost::lexical_cast<std::string>(r.k())
                    + ".");
            }
            if (prob > rnd)
            {
                LOG_DEBUG(("fire reaction"));
                const typename reaction_rule_type::species_id_range products(
                    r.get_products());

                switch (::size(products))
                {
                case 1:
                    {
                        const species_id_type product(products[0]);
                        const species_type sp(tx_.get_species(product));

                        const position_type new_pos(
                            tx_.apply_boundary(
                                divide(
                                    add(multiply(pp0.second.position(), s1.D()),
                                        multiply(tx_.cyclic_transpose(
                                            pp1.second.position(),
                                            pp0.second.position()), s0.D())),
                                    (s0.D() + s1.D()))));
                        boost::scoped_ptr<particle_id_pair_and_distance_list> overlapped(
                            tx_.check_overlap(particle_shape_type(new_pos, sp.radius()),
                                              pp0.first, pp1.first));
                        if (overlapped && overlapped->size() > 0)
                        {
                            throw propagation_error("no space");
                        }

                        if (vc_)
                        {
                            if (!(*vc_)(
                                    particle_shape_type(new_pos, sp.radius()), 
                                    pp0.first, pp1.first))
                            {
                                throw propagation_error("no space");
                            }
                        }

                        remove_particle(pp0.first);
                        remove_particle(pp1.first);
                        particle_id_pair npp(tx_.new_particle(product, new_pos));
                        if (rrec_)
                        {
                            (*rrec_)(
                                reaction_record_type(
                                    r.id(), array_gen(npp.first), pp0.first, pp1.first));
                        }
                        break;
                    }
                case 0:
                    {
                        remove_particle(pp0.first);
                        remove_particle(pp1.first);
                        break;
                    }
                
                default:
                    throw not_implemented("bimolecular reactions that produce more than one product are not supported");
                }

                return true;
            }
        }
        return false;
    }
    
    /*
    bool attempt_reaction(particle_id_pair const& pp, structure_type const& structure)
    {
        reaction_rules const& rules(rules_.query_reaction_rule(pp0.second.sid(), structure_id));
        if (::size(rules) == 0)
        {
            return false;
        }

        const species_type s0(tx_.get_species(pp0.second.sid()));

        const Real rnd(rng_());
        Real prob = 0;

        for (typename boost::range_const_iterator<reaction_rules>::type
                i(boost::begin(rules)), e(boost::end(rules)); i != e; ++i)
        {
            reaction_rule_type const& r(*i);

            boost::shared_ptr<structure_type> pp_structure( 
               tx_.get_structure( tx_.get_species( pp0.second.sid() ).structure_id()) );  

            const Real p( r.k() * dt_ / pp_structure->reaction_volume( s0.radius(), s1.radius(), reaction_length_ ) );
            BOOST_ASSERT(p >= 0.);
            prob += p;
            if (prob >= 1.)
            {
                throw propagation_error(
                    "invalid acceptance ratio ("
                    + boost::lexical_cast<std::string>(p)
                    + ") for reaction rate "
                    + boost::lexical_cast<std::string>(r.k())
                    + ".");
            }
            if (prob > rnd)
            {
                LOG_DEBUG(("fire reaction"));
                const typename reaction_rule_type::species_id_range products(
                    r.get_products());

                switch (::size(products))
                {
                case 1:
                    {
                        const species_id_type product(products[0]);
                        const species_type sp(tx_.get_species(product));

                        const position_type new_pos(
                            tx_.apply_boundary(
                                divide(
                                    add(multiply(pp0.second.position(), s1.D()),
                                        multiply(tx_.cyclic_transpose(
                                            pp1.second.position(),
                                            pp0.second.position()), s0.D())),
                                    (s0.D() + s1.D()))));
                        boost::scoped_ptr<particle_id_pair_and_distance_list> overlapped(
                            tx_.check_overlap(particle_shape_type(new_pos, sp.radius()),
                                              pp0.first, pp1.first));
                        if (overlapped && overlapped->size() > 0)
                        {
                            throw propagation_error("no space");
                        }

                        if (vc_)
                        {
                            if (!(*vc_)(
                                    particle_shape_type(new_pos, sp.radius()), 
                                    pp0.first, pp1.first))
                            {
                                throw propagation_error("no space");
                            }
                        }

                        remove_particle(pp.first);
                        particle_id_pair npp(tx_.new_particle(product, new_pos));
                        if (rrec_)
                        {
                            (*rrec_)(
                                reaction_record_type(
                                    r.id(), array_gen(npp.first), pp0.first, pp1.first));
                        }
                        break;
                    }
                case 0:
                    {
                        remove_particle(pp.first);
                        break;
                    }
                
                default:
                    throw not_implemented("surface interactions that produce more than one product are not supported");
                }

                return true;
            }
        }
        return false;
    }
    */
    
    void remove_particle(particle_id_type const& pid)
    {
        LOG_DEBUG(("remove particle %s", boost::lexical_cast<std::string>(pid).c_str()));
        tx_.remove_particle(pid);
        typename particle_id_vector_type::iterator i(
            std::find(queue_.begin(), queue_.end(), pid));
        if (queue_.end() != i)
            queue_.erase(i);
    }

private:
    position_type random_unit_vector()
    {
        position_type v(rng_() - .5, rng_() - .5, rng_() - .5);
        return v / length(v);
    }

private:
    particle_container_type& tx_;
    network_rules_type const& rules_;
    rng_type& rng_;
    Real const dt_;
    int const max_retry_count_;
    reaction_recorder_type* const rrec_;
    volume_clearer_type* const vc_;
    particle_id_vector_type queue_;
    int rejected_move_count_;
    Real const reaction_length_;
    static Logger& log_;
};

template<typename Ttraits_>
Logger& newBDPropagator<Ttraits_>::log_(Logger::get_logger("ecell.newBDPropagator"));

#endif /* NEW_BD_PROPAGATOR_HPP */

