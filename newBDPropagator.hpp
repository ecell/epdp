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
    typedef std::pair<position_type, position_type> position_pair_type;
    typedef std::pair<position_type, length_type> projected_type;
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

        //LOG_DEBUG(("propagating particle %s", boost::lexical_cast<std::string>(pp.first).c_str()));



        /* START SOME SAMPLES, after which we return. 
        position_pair_type pp01;
        position_pair_type new_pp01;
        
        species_type const sp(tx_.get_species(pp.second.sid()));
        boost::shared_ptr<structure_type> const pp_struct( tx_.get_structure( sp.structure_id() ) );
        
        reaction_rules const& rules(rules_.query_reaction_rule(pp.second.sid()));

        typename boost::range_const_iterator<reaction_rules>::type i(boost::begin(rules));
        reaction_rule_type const& r(*i);
        
        typename reaction_rule_type::species_id_range products(r.get_products());
        
        species_type s0(tx_.get_species(products[0])), s1(tx_.get_species(products[1]));
              
        length_type r01( s0.radius() + s1.radius() );
              
        if(s0.structure_id() == "world")    
            std::swap(s0,s1);

        for(int i = 0; i < 2000; i++)
        {
            pp01 = pp_struct->
                special_geminate_dissociation_positions( rng_, s0, s1, pp.second.position(), reaction_length_ ); 
            
            pp01.first = tx_.apply_boundary( pp01.first );
            pp01.second = tx_.apply_boundary( pp01.second );
            
            if( tx_.distance(pp01.first, pp01.second) - r01 < 0 )
                std::cerr << "Diss pp dist < 0 error" << std::endl;
                
            if( tx_.distance( pp_struct->projected_point( pp01.second ).first, pp01.second) - s1.radius() < 0 )
                std::cerr << "Diss sp dist < 0 error" << std::endl;
                                
            new_pp01 = make_pair_move(s0, s1, pp01, pp.first);           
             
            std::cout << new_pp01.first[0] << " " << new_pp01.first[1] << " " << new_pp01.first[2] << " " 
                      << new_pp01.second[0] << " " << new_pp01.second[1] << " " << new_pp01.second[2] << std::endl;
                      
            if(tx_.distance(new_pp01.first, new_pp01.second) - r01 < 0)
                std::cerr << "Move pp dist < 0 error" << std::endl;
                
            if(tx_.distance(pp_struct->projected_point( new_pp01.second ).first, new_pp01.second) - s1.radius() < 0)
                std::cerr << "Move sp dist < 0 error" << std::endl;
        }     
                  
        return false;  
        */
        
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
            position_type const displacement( pp_structure->
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
        int j = 0;
        while(!bounced && j < particles_in_overlap)
            bounced = overlapped->at(j++).second < r0;
        
        /* If the particle has not bounced with another particle and lives in the bulk: 
           check for core overlap with a surface. */
        if(!bounced && pp_structure->id() == "world")
        {
            /* Check for overlap with nearest surface. */
            struct_id_distance_pair = tx_.get_closest_surface( new_pos );
            core_surface_distance = struct_id_distance_pair.second - r0;
        
            if(core_surface_distance < 0)
                bounced = true;
        }         
         
        /* If particle is bounced, check reaction_volume for reaction partners at old_pos. */     
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
        
        /* Attempt a reaction or interaction with all the particles (and/or a surface) that 
           are inside the reaction volume */
        Real accumulated_prob = 0;
        
        /* First, if a surface is near, attempt an interaction.  */
        if(core_surface_distance < reaction_length_ && pp_structure->id() == "world")
        {
            try
            {
                const Real prob( attempt_interaction(pp, tx_.get_structure(struct_id_distance_pair.first) ) );
                
                if(prob < 0)
                {
                    //An interaction occured.
                    return true;                                 
                }
                else                
                {
                    LOG_DEBUG(("Particle attempted an interaction with the noninteractive surface %s.", 
                        boost::lexical_cast<std::string>(struct_id_distance_pair.first).c_str()));
                    accumulated_prob += prob;
                    ++rejected_move_count_;
                }
            }
            catch (propagation_error const& reason)
            {
                log_.info("surface interaction rejected (reason: %s)", reason.what());
                ++rejected_move_count_;
            }
        }
        
        /* Now attempt reaction with all particles inside the reaction volume. */
        j = 0;
        while(j < particles_in_overlap)
        {
            particle_id_pair_and_distance const& closest( overlapped->at(j) );
 		        
    		try
            {
                const Real prob( attempt_reaction(pp, closest.first) );
            
                if( prob < 0 )
                    //A reaction occured.
                    return true;
                else                
                {
                    LOG_DEBUG(("Particle attempted a reaction with a nonreactive particle %s.", 
                        boost::lexical_cast<std::string>(closest.first.first).c_str()));
                    accumulated_prob += prob;
                    ++rejected_move_count_;
                }
            }
            catch (propagation_error const& reason)
            {
                log_.info("second-order reaction rejected (reason: %s)", reason.what());
                ++rejected_move_count_;
            }
            
            j++;
        }
                       
        /* If the particle did neither react, interact or bounce, update it to it's new position. */
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
                        boost::shared_ptr<structure_type> pp_struct( tx_.get_structure( tx_.get_species(pp.second.sid()).structure_id() ) );
                        position_type new_pos = pp.second.position();
                        
                        /* If the product particle does NOT live of the same structure as pp => surface -> bulk dissociation. */
                        if( s0.structure_id() != pp_struct->id() )
                        {
                            if( s0.structure_id() != "world" )
                                throw not_implemented("No surface -> surface dissociation allowed");
                        
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
                            throw propagation_error("no space due to particle");
                        }

                        if (vc_)
                        {
                            if (!(*vc_)(new_p.second.shape(), pp.first))
                            {
                                throw propagation_error("no space due to particle");
                            }
                        }

                        if( s0.structure_id() == "world" )
                        {
                            structure_id_and_distance_pair const struct_id_distance_pair( tx_.get_closest_surface( new_pos ) );
                            length_type const core_surface_distance( struct_id_distance_pair.second - s0.radius() );
                                    
                            if( core_surface_distance < 0 )
                                throw propagation_error("no space due to surface");
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
                        species_type s0(tx_.get_species(products[0])),
                                s1(tx_.get_species(products[1]));
                        int i = max_retry_count_;
                        position_pair_type pp01;

                        boost::shared_ptr<structure_type> pp_struct( tx_.get_structure( tx_.get_species( pp.second.sid() ).structure_id() ) );
                        
                        /* Create positions (np0, np1) for the reaction products. 
                        The ipv between the products lies within the reation volume. */
                        for (;;)
                        {
                            if (--i < 0)
                            {
                                throw propagation_error("no space");
                            }
                            
                            /* If the product particles do NOT live on the same structure => surface -> surface + bulk dissociation.
                               Else, the products live on the same structure AND on the same structure as the reactant. 
                               => standard geminate dissociation. 
                               
                               We set: species 0 lives on the surface; species 1 lives in the bulk. */
                            if( s0.structure_id() != s1.structure_id() )
                            {
                                if( !(s0.structure_id() == "world" || s1.structure_id() == "world")  )
                                    throw not_implemented("No surface -> surface + surface - dissociation allowed");
                                                 
                                if(s0.structure_id() == "world")    
                                     std::swap(s0,s1);

                                pp01 = pp_struct->
                                        special_geminate_dissociation_positions( rng_, s0, s1, pp.second.position(), reaction_length_ );
                            }
                            else
                            {
                                pp01 = pp_struct->
                                        geminate_dissociation_positions( rng_, s0, s1, pp.second.position(), reaction_length_ );
                            }
                            
                            pp01.first = tx_.apply_boundary( pp01.first );
                            pp01.second = tx_.apply_boundary( pp01.second );
                                                                                   
                            boost::scoped_ptr<particle_id_pair_and_distance_list> overlapped_s0(
                                tx_.check_overlap(
                                    particle_shape_type(pp01.first, s0.radius()),
                                    pp.first));
                            boost::scoped_ptr<particle_id_pair_and_distance_list> overlapped_s1(
                                tx_.check_overlap(
                                    particle_shape_type(pp01.second, s1.radius()),
                                    pp.first));
                                    
                            //Only when both products live in the bulk, check for possible overlap with a surface.
                            bool surface_overlap = false;
                            if( s0.structure_id() == "world" && s1.structure_id() == "world" )
                            {
                                structure_id_and_distance_pair struct_id_distance_pair( tx_.get_closest_surface( pp01.first ) );
                                length_type core_surface_distance( struct_id_distance_pair.second - s0.radius() );
                                    
                                if( core_surface_distance < 0 )
                                {
                                    surface_overlap = true;
                                }
                                    
                                struct_id_distance_pair = tx_.get_closest_surface( pp01.second );
                                core_surface_distance = struct_id_distance_pair.second - s1.radius();

                                if( core_surface_distance < 0 )
                                    surface_overlap = true;                                
                            }
                                    
                            if (!(overlapped_s0 && overlapped_s0->size() > 0) && !(overlapped_s1 && overlapped_s1->size() > 0) && !surface_overlap)
                                break;
                        }

                        /* After dissociation both particles are allowed to move. */
                        pp01 = make_pair_move(s0, s1, pp01, pp.first);

                        if (vc_)
                        {
                            if (!(*vc_)(particle_shape_type(pp01.first, s0.radius()), pp.first) || !(*vc_)(particle_shape_type(pp01.second, s1.radius()), pp.first))
                            {
                                throw propagation_error("no space");
                            }
                        }

                        tx_.remove_particle(pp.first);
                        const particle_id_pair
                            npp0(tx_.new_particle(s0.id(), pp01.first)),
                            npp1(tx_.new_particle(s1.id(), pp01.second));

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

    /* Function returns the accumulated acceptance probability of the reacting particles if no
       reaction occured. If a reaction has happened, it returns -1 and if no rules exist it returns 0. */
    Real attempt_reaction(particle_id_pair const& pp0, particle_id_pair const& pp1)
    {
        reaction_rules const& rules(rules_.query_reaction_rule(pp0.second.sid(), pp1.second.sid()));
        if (::size(rules) == 0)
        {
            return 0;
        }

        species_type s0(tx_.get_species(pp0.second.sid())),
                s1(tx_.get_species(pp1.second.sid()));
        
        /* If the structures of the reactants do not equal, one of the reactants has to come from the bulk,
           and we let this be s1, the particle from the surface is named s0. */     
        if(s0.structure_id() != s1.structure_id())
        {
            if( !(s0.structure_id() == "world" || s1.structure_id() == "world")  )
                throw not_implemented("No surface -> surface + surface - reactions allowed");
                
            if(s0.structure_id() == "world")
                std::swap(s0,s1);
        }

        const Real rnd(rng_());
        Real prob = 0;

        boost::shared_ptr<structure_type> s1_struct( tx_.get_structure( s1.structure_id()) );
        const Real reaction_volume( s1_struct->particle_reaction_volume( s0.radius() + s1.radius(), reaction_length_ ) );
        
        for (typename boost::range_const_iterator<reaction_rules>::type
                i(boost::begin(rules)), e(boost::end(rules)); i != e; ++i)
        {
            reaction_rule_type const& r(*i);  

            //TODO Do I need special reaction volumes for surface-bulk interaction? Or allways the 3D case?
            /* NOTE: multiplying reaction volume with 2 because a reaction between a pair of particles is
               attempted twice (at max). This only holds for small Pacc. */
            const Real p( r.k() * dt_ / ( 2 * reaction_volume ) );
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
                        
                        position_type new_pos( tx_.apply_boundary(
                                        divide(
                                        add(multiply(pp0.second.position(), s1.D()),
                                            multiply(tx_.cyclic_transpose(
                                                pp1.second.position(),
                                                pp0.second.position()), s0.D())),
                                        (s0.D() + s1.D()))) );

                        //For unequal structures, project new_pos on surface.
                        if(sp.structure_id() != s1.structure_id())
                        {
                            boost::shared_ptr<structure_type> sp_struct( tx_.get_structure( sp.structure_id() ) );
                            projected_type const position_on_surface( sp_struct->projected_point(new_pos) );
                            
                            new_pos = tx_.apply_boundary( position_on_surface.first );
                        }

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
                        
                        if( sp.structure_id() == "world" )
                        {
                            structure_id_and_distance_pair const struct_id_distance_pair( tx_.get_closest_surface( new_pos ) );
                            length_type const core_surface_distance( struct_id_distance_pair.second - sp.radius() );
                                    
                            if( core_surface_distance < 0 )
                                throw propagation_error("no space due to surface");     
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
                
                return -1;
            }
        }
        return prob;
    }
    
    /* Function returns the accumulated acceptance probability of the reacting particles if no
       reaction occured. If a reaction has happened, it returns -1 and if no rules exist it returns 0. */
    Real attempt_interaction(particle_id_pair const& pp, boost::shared_ptr<structure_type> const& structure)
    {
        reaction_rules const& rules(rules_.query_reaction_rule(pp.second.sid(), structure->sid() ));
        if (::size(rules) == 0)
        {
            return 0;
        }
        
        const species_type s0(tx_.get_species(pp.second.sid()));
        
        boost::shared_ptr<structure_type> pp_structure( tx_.get_structure( s0.structure_id() ) ); 

        if(pp_structure->id() != "world")
            throw not_implemented("only particles living in the bulk are allowed to interact with surfaces");
        
        const Real rnd(rng_());
        Real prob = 0;

        for (typename boost::range_const_iterator<reaction_rules>::type
                i(boost::begin(rules)), e(boost::end(rules)); i != e; ++i)
        {
            reaction_rule_type const& r(*i); 

            const Real p( r.k() * dt_ / structure->surface_reaction_volume( s0.radius(), reaction_length_ ) );
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
                LOG_DEBUG(("fire surface interaction"));
                const typename reaction_rule_type::species_id_range products(
                    r.get_products());

                switch (::size(products))
                {
                case 1:
                    {
                        const species_id_type product(products[0]);
                        const species_type sp(tx_.get_species(product));

                        const position_type new_pos( tx_.apply_boundary( structure->projected_point( pp.second.position() ).first ) );
                            
                        boost::scoped_ptr<particle_id_pair_and_distance_list> overlapped(
                            tx_.check_overlap(particle_shape_type(new_pos, sp.radius()),
                                              pp.first));
                        if (overlapped && overlapped->size() > 0)
                        {
                            throw propagation_error("no space");
                        }

                        if (vc_)
                        {
                            if (!(*vc_)(
                                    particle_shape_type(new_pos, sp.radius()), 
                                    pp.first))
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
                                    r.id(), array_gen(npp.first), pp.first));
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

                return -1;
            }
        }
        return prob;
    }
    
    /* Given a position pair and two species, the function generates two new 
       positions based on two moves drawn from the free propegator. */
    position_pair_type const make_pair_move(species_type const& s0, species_type const& s1, position_pair_type const& pp01, 
                                                particle_id_type const& ignore)
    {
        boost::shared_ptr<structure_type> s0_struct( tx_.get_structure( s0.structure_id() ) );
        boost::shared_ptr<structure_type> s1_struct( tx_.get_structure( s1.structure_id() ) );

        length_type const r01( s0.radius() + s1.radius() );
        
        position_type temp_pos0, temp_pos1;
        position_pair_type new_pp01( pp01 );
        
        //Randomize which species moves first, and if there will be one or two moves. (50%)
        bool move_prt0_first( rng_.uniform_int(0, 1) );
        int move_count( 1 + rng_.uniform_int(0, 1) );
        while( move_count-- )
        {
            if( move_prt0_first )
            {
                temp_pos0 = make_move(s0, pp01.first, ignore);
                new_pp01.first = tx_.distance(temp_pos0, new_pp01.second) - r01 < 0 ? pp01.first : temp_pos0;
            }                       
            else
            {
                temp_pos1 = make_move(s1, pp01.second, ignore);
                new_pp01.second = tx_.distance(temp_pos1, new_pp01.first) - r01 < 0 ? pp01.second : temp_pos1;
            }
            move_prt0_first = !move_prt0_first;
        }  
            
        return new_pp01;
    }

    position_type const make_move(species_type const& s0, position_type const& np, particle_id_type const& ignore)
    {
        if(s0.D() == 0)
            return np;
    
        boost::shared_ptr<structure_type> s0_struct( tx_.get_structure( s0.structure_id() ) );

        position_type displacement( s0_struct->
                                        bd_displacement(s0.v() * dt_, std::sqrt(2.0 * s0.D() * dt_), rng_) );

        position_type new_pos(tx_.apply_boundary( add( np, displacement ) ));
                        
        boost::scoped_ptr<particle_id_pair_and_distance_list> overlapped( 
                                    tx_.check_overlap(
                                    particle_shape_type(new_pos, s0.radius()),
                                    ignore));
                                                                
        if ( (overlapped && overlapped->size() > 0) )
            return np;
        
        if(s0.structure_id() == "world")
        {
            structure_id_and_distance_pair const struct_id_distance_pair( tx_.get_closest_surface( new_pos ) );
            length_type const core_surface_distance( struct_id_distance_pair.second - s0.radius() );
                                    
            if( core_surface_distance < 0 )
                return np;        
        }
            
        return new_pos;
    }
    
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

