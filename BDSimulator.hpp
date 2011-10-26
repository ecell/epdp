#ifndef BD_SIMULATOR_HPP
#define BD_SIMULATOR_HPP

#include <algorithm>
#include <limits>
#include <boost/foreach.hpp>
#include "NetworkRules.hpp"
#include "newBDPropagator.hpp"
#include "World.hpp"
#include "ParticleSimulator.hpp"
#include "utils/pair.hpp"

template<typename Tworld_>
struct BDSimulatorTraitsBase: public ParticleSimulatorTraitsBase<Tworld_>
{
};

template<typename Ttraits_>
class BDSimulator: public ParticleSimulator<Ttraits_>
{
public:
    typedef Ttraits_ traits_type;
    typedef ParticleSimulator<Ttraits_> base_type;
    typedef typename traits_type::world_type world_type;
    typedef typename world_type::traits_type::rng_type rng_type;
    typedef typename world_type::species_id_type species_id_type;
    typedef typename world_type::species_type species_type;
    typedef typename world_type::length_type length_type;
    typedef typename traits_type::time_type time_type;  
    typedef typename traits_type::network_rules_type network_rules_type;
    typedef typename traits_type::reaction_rule_type reaction_rule_type;
    typedef typename traits_type::rate_type rate_type;
    typedef typename traits_type::reaction_record_type reaction_record_type;
    typedef typename traits_type::reaction_recorder_type reaction_recorder_type;

public:
    Real const& dt_factor()
    {
        return dt_factor_;
    }

    virtual ~BDSimulator() {}

    BDSimulator(boost::shared_ptr<world_type> world, 
                boost::shared_ptr<network_rules_type const> network_rules,
                rng_type& rng, Real dt_factor = traits_type::DEFAULT_DT_FACTOR,
                int dissociation_retry_moves = 1, length_type reaction_length_factor = 1)
        : base_type(world, network_rules, rng),
          dt_factor_(dt_factor), num_retries_(dissociation_retry_moves), 
          reaction_length_factor_(reaction_length_factor)
    {
        calculate_dt();
        reaction_length_ = determine_reaction_length( *world, base_type::dt_, reaction_length_factor_ );
        bound_time = 0;
    }

    virtual void calculate_dt()
    {
        base_type::dt_ = dt_factor_ * determine_dt(*base_type::world_);
        LOG_DEBUG(("dt=%f", base_type::dt_));
    }

    virtual void step()
    {
        _step(base_type::dt_);
    }

    virtual bool step(time_type upto)
    {
        time_type const lt(upto - base_type::t_);
        if (lt <= 0.)
            return false;
        if (base_type::dt_ < lt)
        {
            _step(base_type::dt_);
        }
        else
        {
            _step(lt);
            base_type::t_ = upto;
        }
        
        if(base_type::world_->num_particles() == 1)
            bound_time += base_type::dt_;
        
        return true;
    }

    static Real determine_dt(world_type const& world)
    {
        Real D_max(0.), radius_min(std::numeric_limits<Real>::max());

        BOOST_FOREACH(species_type s, world.get_species())
        {
            if (D_max < s.D())
                D_max = s.D();
            if (radius_min > s.radius())
                radius_min = s.radius();
        }
        return gsl_pow_2(radius_min * 2) / (D_max * 2);
    }
    
    static length_type determine_reaction_length(world_type const& world, time_type dt, length_type const& reaction_length_factor)
    {
        return reaction_length_factor * sqrt( 2 * 1E-12 * dt );
    }
    
    Real get_reaction_length() const
    {
        return reaction_length_;    
    }
    
    double set_reaction_length_factor(length_type new_reaction_length_factor)
    {
        double old_bound_time = bound_time;
        bound_time = 0;
        
        reaction_length_factor_ = new_reaction_length_factor;
        reaction_length_ = determine_reaction_length( *base_type::world_, base_type::dt_, reaction_length_factor_ );
        return old_bound_time;        
    }

protected:
    void _step(time_type dt)
    {
        {
            newBDPropagator<traits_type> propagator(
                *base_type::world_,
                *base_type::network_rules_,
                base_type::rng_,
                dt, num_retries_, reaction_length_,
                base_type::rrec_.get(), 0,
                make_select_first_range(base_type::world_->
                                        get_particles_range()));
            while (propagator());
            //LOG_DEBUG(("%d: t=%lg, dt=%lg", base_type::num_steps_, 
            //           base_type::t_, dt));
        }
        ++base_type::num_steps_;
        base_type::t_ += dt;
    }

private:
    Real const dt_factor_;
    int const num_retries_;
    length_type reaction_length_factor_;
    length_type reaction_length_;
    static Logger& log_;
    Real bound_time;
};

template<typename Ttraits_>
Logger& BDSimulator<Ttraits_>::log_(Logger::get_logger("BDSimulator"));


#endif /* BD_SIMULATOR_HPP */
