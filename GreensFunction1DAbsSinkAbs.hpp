#if !defined( __GREENSFUNCTION1DABSSINKABS_HPP )
#define __GREENSFUNCTION1DABSSINKABS_HPP

#include <sstream>
#include <iostream>
#include <cstdlib>
#include <exception>
#include <vector>

#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_trig.h>
#include <gsl/gsl_sum.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_sf_expint.h>
#include <gsl/gsl_sf_elljac.h>
#include <gsl/gsl_roots.h>

#include <math.h>

#include "findRoot.hpp"
#include "Defs.hpp"
#include "OldDefs.hpp"			// TODO: this must be removed at some point!
#include "GreensFunction.hpp"

class GreensFunction1DAbsSinkAbs: public GreensFunction
{
private:
    // This is a typical length scale of the system, may not be true!
    static const Real L_TYPICAL = 1E-8;
    // The typical timescale of the system, may also not be true!!
    static const Real T_TYPICAL = 1E-6;
    // measure of 'sameness' when comparing floating points numbers
    static const Real EPSILON = 1E-10;
    // Is 1E3 a good measure for the probability density?!
    static const Real PDENS_TYPICAL = 1;
    // The maximum number of terms used in calculating the sum
    static const int MAX_TERMS = 500;
    // The minimum number of terms
    static const int MIN_TERMS = 20;

public:
    GreensFunction1DAbsSinkAbs(Real D, Real k, Real r0, Real rsink, Real sigma, Real a)
	: GreensFunction(D), v(0.0), k(k), r0(r0), sigma(sigma), a(a), rsink(rsink), l_scale(L_TYPICAL), t_scale(T_TYPICAL)
    {
	    /* Set variables which define a domain with the sink at the origin. 
	       Futhermore r0 is assumed to be ringht from the sink. */
	    assert( a > sigma );
	    
	    L0 = fabs( r0 - rsink ); 
	    if( r0 >= rsink )
	    {
	        Lr = a - rsink; 
	        Ll = rsink - sigma;	        
	    }
	    else
	    {
	        Lr = rsink - sigma; 
	        Ll = a - rsink;
	    }
    }

    /* The constructor is overloaded and can be called with or without drift v
       copy constructor including drift variable v */
    GreensFunction1DAbsSinkAbs(Real D, Real v, Real k, Real r0, Real rsink, Real sigma, Real a)
	: GreensFunction(D), v(v), k(k), r0(r0), sigma(sigma), a(a), rsink(a - sigma), l_scale(L_TYPICAL), t_scale(T_TYPICAL)
    {
		assert( a > sigma );
	    
	    L0 = fabs( r0 - rsink ); 
	    if( r0 >= rsink )
	    {
	        Lr = a - rsink; 
	        Ll = rsink - sigma;	        
	    }
	    else
	    {
	        Lr = rsink - sigma; 
	        Ll = a - rsink;
	    }
    }

    ~GreensFunction1DAbsSinkAbs()
    {
	;   // empty
    }

    Real geta() const
    {
	    return a;
    }
    
    Real getLr() const
    {
	    return Lr;
    }
    
    Real getLl() const
    {
	    return Ll;
    }
    
    Real getL0() const
    {
	    return L0;
    }

    Real getrsink() const
    {
	    return rsink;
    }
    
    Real getsigma() const
    {
	    return sigma;
    }

    Real getr0() const
    {
	    return r0;
    }

    Real getk() const
    {
	    return k;
    }

    Real getv() const
    {
	    return v;
    }

    
    /* Calculates the probability density of finding the particle at 
       location z at timepoint t, given that the particle is still in the 
       domain. */
    Real calcpcum(Real r, Real t) const;

    /* Determine which event has occured at time t. Either an escape 
       (left or right boundary) or a reaction with the sink. 
       Based on the fluxes through the boundaries at the given time. */
    EventKind drawEventType( Real rnd, Real t ) const;

    /* Draws the first passage time from the propensity function */
    Real drawTime(Real rnd) const;

    /* Draws the position of the particle at a given time, assuming that 
       the particle is still in the domain. */
    Real drawR(Real rnd, Real t) const;

    /* Calculates the probability flux leaving the domain through the right 
       absorbing boundary at time t. */
    Real flux_leavea(Real t) const;

    /* Calculates the probability flux leaving the domain through the left 
       absorbing boundary at time t. */
    Real flux_leaves(Real t) const;


    /* TODO: Private methods which are public for now - debugging */

    // Calculates the probability of finding the particle inside the 
    // domain at time t -> the survival probability
    Real p_survival(Real t) const;

    // Calculates the total probability flux leaving the domain at time t
    Real flux_tot(Real t) const;

    /* Calculates the probability flux leaving the domain through the 
       sink, sub-domain containing r0 via absorbing boundary and sub-domain
       not containing r0 via sink, respectivily. */
    Real flux_sink(Real t) const;
    Real flux_abs_Lr(Real t) const;
    Real flux_abs_Ll(Real t) const;

    // Calculates the probability density of finding the particle at 
    // location r at time t.
    Real prob_r(Real r, Real t) const;
    
// End of public/private mix methods

//private:	// method made public for testing

    std::string dump() const;

    const char* getName() const
    {
        return "GreensFunction1DAbsSinkAbs";
    }

    // Calculates the roots of x * sin(x) + k * L / (2 * D) * ( cos(x * (Lr - Ll)/L) - cos(x) ) 
    Real root_n(int const& n) const;
    
private:

    struct drawR_params
    {
	    Real root_n[MAX_TERMS];
	    Real exp_and_denominator[MAX_TERMS];
	    Real prefactor;
        Real L0;
        Real Lr;
        Real Ll;
	    Real rnd;
        Real D;
        Real k;
	    int terms;
    };

    struct drawT_params
    {
	    Real exponent[MAX_TERMS];
	    Real Xn[MAX_TERMS];
	    Real prefactor;
	    int  terms;
	    // the timescale used for convergence
	    Real   tscale;
	    Real rnd;
    };
    
    struct root_f_params
    {
	    Real Lm_L;
	    Real h;
    };
       
    inline Real p_denominator( Real const& root_n ) const;
    
    inline Real Greens_fn(Real const& t, Real const& root_n, Real const& root_n2) const;

    static Real num_int_r_leftdomain(Real const& rr, Real const& root_n, drawR_params const& params);

    static Real num_int_r_rightdomainA(Real const& rr, Real const& root_n, drawR_params const& params);

    static Real num_int_r_rightdomainB(Real const& rr, Real const& root_n, drawR_params const& params);

    static Real root_f (Real x, void *p);

    static Real drawT_f (Real t, void *p);

    static Real drawR_f (Real t, void *p);

    /* Class variables */
    // The diffusion constant and drift velocity
    const Real v;
    // The reaction constant
    const Real k;
    //starting position
    const Real r0;
    // The left and right boundary of the domain (sets the l_scale, see below)
    const Real sigma;
    const Real a;
    //Position of the sink in the domain.
    const Real rsink;
    // This is the length scale of the system
    Real l_scale;
    // This is the time scale of the system.
    Real t_scale;
    
    /* Greensfunction assumes that the sink is at the origin, and
       consists of two sub-domains: one between a boundary and the sink including
       r0, and one between boundary and sink not inlcuding r0. */

    // Length of sub-domain which does not include r0.
    Real Lr;
    // Length of sub-domain which does include r0.
    Real Ll;
    // Distance between the sink and r0.
    Real L0;
};
#endif // __GREENSFUNCTION1DRADABS_HPP
