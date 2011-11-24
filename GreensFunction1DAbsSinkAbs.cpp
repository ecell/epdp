#include <sstream>
#include <iostream>
#include <cstdlib>
#include <exception>
#include <vector>
#include <math.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_trig.h>
#include <gsl/gsl_sum.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_sf_expint.h>
#include <gsl/gsl_sf_elljac.h>
#include <gsl/gsl_roots.h>

#include "findRoot.hpp"
#include "GreensFunction1DAbsSinkAbs.hpp"

/* This is the appropriate definition of the function defining
   the roots of our Green's functions in GSL.
   Later needed by the rootfinder. */
Real GreensFunction1DAbsSinkAbs::root_f (Real x, void *p)
{
    struct root_f_params *params = (struct root_f_params *)p;
    const Real Lm_L = (params->Lm_L);
    const Real h = (params->h);

    // L   = Lr + Ll
    // h    = L * k / (2 * D)
    // L_Lm = Lr + Ll / Lr - Ll
    // x    = q * L

    return x * sin(x) + h * ( cos(x * Lm_L) - cos(x) );

}

/* Calculates the first n roots of root_f */
void GreensFunction1DAbsSinkAbs::calculate_n_roots(int const& n) const
{
    const Real Lr( getLr() );
    const Real Ll( getLl() );
    const Real L( Lr + Lr );
    const Real Lm( Lr - Ll );
    const Real Lm_L( Lm / L );
    const Real h( getk() * L / ( 2 * getD() ) );
    
    Real root_i;
    unsigned int i = 0;

    /* Stores the last roots found with a low and a high
       frequency resp. */
    real_pair last_roots, lower_upper_pair;

    /* set params needed for get lower and upper function */
    lower_upper_params lo_up_params;
    lo_up_params.h = h;
    lo_up_params.Lm_L = Lm_L;
    lo_up_params.long_half_wave = max( L/Lr * M_PI, L/Ll * M_PI );
    lo_up_params.short_half_wave = min( L/Lr * M_PI, L/Ll * M_PI );

    /* Define the root function. */
    gsl_function F;
    struct root_f_params params = { Lm_L, h };
     
    F.function = &GreensFunction1DAbsSinkAbs::root_f;
    F.params = &params;

    /* define and ad a new solver type brent */
    const gsl_root_fsolver_type* solverType( gsl_root_fsolver_brent );
    gsl_root_fsolver* solver( gsl_root_fsolver_alloc( solverType ) );

    /* Find the n roots */
    while(i++ < n)
    {
        lower_upper_pair = get_lower_and_upper( last_roots, lower_upper_params );

        root_i = findRoot( F, solver, lower_upper_pair.first, lower_upper_pair.second, 
                           1.0*EPSILON, EPSILON, "GreensFunction1DAbsSinkAbs::root_f" );

        assert( root_i > get_last_root() - EPSILON );

        ad_to_rootList( root_i );       

        if( root_i - last_roots.first < root_i - last_roots.second )
            last_roots.first = root_i;
        else
            last_roots.second = root_i;
    }

    gsl_root_fsolver_free( solver );
}

/* Function returns two points on the x-axis which straddle the next root. */
real_pair GreensFunction1DAbsSinkAbs::get_lower_and_upper( real_pair& last_roots, 
                                                           lower_upper_params const& params )
{
    const Real root_n( get_last_root() );
    const Real safety( .75 );

    Real lower, upper, next_root_est, max_offset;

    if( params.h / root_n < 1 )
    {
        max_offset = M_PI;
        next_root_est = root_n + M_PI;
    }
    else
    {
        const Real next_root_long( last_roots.first + params.long_half_wave );
        const Real next_root_short( last_roots.second + params.short_half_wave );

        if( next_root_long < next_root_short )
        {
            next_root_est = next_root_long;

            max_offset = min( next_root_short - next_root_est, 
                              params.long_half_wave );
        }
        else
        {
            next_root_est = root_n + next_root_short;
            
            max_offset = min( next_root_long - next_root_est, 
                              params.short_half_wave );
        }
    }
       
    lower = next_root_est - safety * max_offset;
    upper = next_root_est + safety * max_offset;

    struct root_f_params p = { params.Lm_L, params.h };

    Real f_lower( root_f( lower, p ) );
    Real f_upper( root_f( upper, p ) );
    bool odd( (get_rootList.size() + 1)%2 );
    /* if 'endpoints do not straddle' */
    if( f_lower * f_higher > 0 )
    {
        std::cerr << "An endpoints do not straddle is imminent!" << std::endl;
        std::cerr << "f_low(" << lower << ") = " << f_lower << ", f_high(" 
                  << upper << ") = " << f_upper << std::endl;
        /*
        if( odd )
        {
            while( f_lower < 0 )
            {
                
            }
                

        }
        else
        {

        }
        */
    }    

    return real_pair( lower, upper );
}

/* returns a guess for the number of terms needed for 
   the greensfunction to converge at time t */
unsigned int GreensFunction1DAbsSinkAbs::guess_maxi(Real t) const
{
    const unsigned int safety(2);

    if (t >= INFINITY)
    {
        return safety;
    }

    const Real D( getD() );

    const Real root0( get_root( 0 ) );
    const Real Dt(D * t);

    const Real thr(exp(- Dt * root0 * root0) * this->TOLERANCE * 1e-1);

    if (thr <= 0.0)
    {
        return this->MAX_TERMS;
    }

    const Real max_root(sqrt(root0 * root0 - log(thr) / Dt));

    const unsigned int 
        maxi(safety + 
             static_cast<unsigned int>(max_alpha * getL() / M_PI));

    return std::min(maxi, MAX_TERMS);
}

/* Standart form of the greensfunction */
inline Real GreensFunction1DAbsSinkAbs::Greens_fn(Real const& t, 
                                                  Real const& root_n, 
                                                  Real const& root_n2) const
{
    return exp( - getD() * root_n2 * t ) / p_denominator( root_n );
}

/* Denominator of the greensfunction. */
inline Real GreensFunction1DAbsSinkAbs::p_denominator(Real const& root_n) const
{
    const Real Lm( getLr() - getLl() );
    const Real L( getLr() + getLl() );
    
    const Real term1( root_n * L * cos( root_n * L ) + sin( root_n * L ) );
    const Real term2( L * sin( root_n * L ) - Lm * sin( root_n * Lm ) );   

    return getD() * term1 + getk() / 2. * term2;
}

/* Calculates the probability of finding the particle inside the domain
   at time t, the survival probability. */
Real GreensFunction1DAbsSinkAbs::p_survival(Real t) const
{
    THROW_UNLESS( std::invalid_argument, t >= 0.0 );
  
    const Real D ( getD() );
    const Real k ( getk() ); 
    const Real Lr( getLr() );
    const Real Ll( getLl() );
    const Real L0( getL0() );
    const Real L ( Lr + Ll );

    if (t == 0.0 || (D == 0.0 && v == 0.0) )
    {
	    //particle can't escape.
	    return 1.0;
    }

    Real sum = 0, numerator = 0, prev_term = 0;
    Real term1, term2, term_n, root_n;
    int n = 1;
    do
    {
        if ( n >= MAX_TERMS )
	    {
	        std::cerr << " Too many terms needed for GF1DSink::p_survival. N: "
	              << n << std::endl;
	        break;
	    }
	    prev_term = numerator;
    
        root_n = get_root(n);
        term1 = sin( root_n * L ) - sin( root_n * (Lr - L0) ) - sin( root_n * (Ll + L0) );
        term2 = sin( root_n * Lr ) - sin( root_n * L0 ) - sin( root_n * (Lr - L0) );
        
        numerator = ( D * term1 + k * sin( root_n * Ll ) * term2 / root_n );
        term_n = Greens_fn(t, root_n, gsl_pow_2( root_n ) ) * numerator;
        //std::cerr << term_n << " " << root_n * L / M_PI << std::endl;
	    sum += term_n;
	    n++;
    }
    while ( fabs(term_n/sum) > EPSILON  ||
	        fabs(prev_term/sum) > EPSILON ||
	        n <= MIN_TERMS );

    return 2.0 * sum;
}


// Calculates the probability density of finding the particle at location r
// at time t.
Real GreensFunction1DAbsSinkAbs::prob_r(Real r, Real t) const
{
    THROW_UNLESS( std::invalid_argument, t >= 0.0 );
    THROW_UNLESS( std::invalid_argument, (r-sigma) >= 0.0 && r <= a && (r0 - sigma) >= 0.0 && r0<=a );
    
    const Real D( getD() );
    const Real k( getk() ); 
    const Real Lr( getLr() );
    const Real Ll( getLl() );
    Real L0( getL0() );
    Real rr( r - getrsink() ) ;

    const Real L( Lr + Ll );

    // if there was no time change or zero diffusivity => no movement
    if (t == 0 || D == 0)
    {
	    // the probability density function is a delta function
	    if (r == r0)
	    {
	        return INFINITY;
	    }
	    else
	    {
	        return 0.0;
	    }
    }

    // if r is at one of the the absorbing boundaries
    if ( fabs( a - r ) < EPSILON * L || fabs( r - sigma ) < EPSILON * L  )
    {
	    return 0.0;
    }

    Real root_n;
    Real sum = 0, term_n = 0, prev_term = 0, numerator = 0;
    int n = 1;
           
    /* Determine wether rr lies in the same sub-domain as r0.
       A different function is caculated when this is the case. */
    if( rr * ( getr0() - getrsink() ) >= 0 )
    {
        if(rr < 0)
            rr *= -1;

        /* Check if r is left or right of the starting position r0.
           If so, interchange r with r0. */
        if( rr < L0 )
        {
            const Real temp( rr );
            rr = L0;
            L0 = temp;
        }

        const Real LlpL0( Ll + L0 );
        const Real Lrmrr( Lr - rr );

        do
        {   
	        if ( n >= MAX_TERMS )
	        {
	            std::cerr << " Too many terms needed for GF1DSink::prob_r. N: "
                          << n << std::endl;
	            break;
	        }
	        prev_term = term_n;

	        root_n = get_root(n);

            numerator =  D * root_n * sin( root_n * LlpL0 ) + 
                k * sin( root_n * Ll ) * sin( root_n * L0);
            numerator *= sin( root_n * Lrmrr );

            term_n = Greens_fn(t, root_n, gsl_pow_2( root_n ) ) * numerator;
            sum += term_n;

    	    n++;
       }
       while ( fabs(term_n/sum) > EPSILON*PDENS_TYPICAL || 
	           fabs(prev_term/sum) > EPSILON*PDENS_TYPICAL ||
	           n <= MIN_TERMS );
	}
    else
    {
        if( rr > 0 )
            rr *= -1;

        const Real LrmL0( Lr - L0 );
        const Real Llprr( Ll + rr );

        do
        {   
	        if ( n >= MAX_TERMS )
	        {
	            std::cerr << " Too many terms needed for GF1DSink::prob_r. N: "
                          << n << std::endl;
	            break;
	        }
	        prev_term = term_n;

	        root_n = this->root_n(n);

            numerator =  D * root_n * sin( root_n * (Llprr) ) * sin( root_n * (LrmL0) ); 
            term_n = Greens_fn(t, root_n, gsl_pow_2( root_n ) ) * numerator;

            sum += term_n;
            
    	    n++;
       }
       while ( fabs(term_n/sum) > EPSILON*PDENS_TYPICAL || 
	           fabs(prev_term/sum) > EPSILON*PDENS_TYPICAL ||
	           n <= MIN_TERMS );
    }
    
    return - 2.0 * sum;
}

/* Calculates the probability density of finding the particle at location z at
   timepoint t, given that the particle is still in the domain. */
Real GreensFunction1DAbsSinkAbs::calcpcum(Real r, Real t) const
{
    return prob_r(r, t)/p_survival(t);
}

/* Function returns flux at absorbing bounday sigma. */
Real GreensFunction1DAbsSinkAbs::flux_leaves(Real t) const
{
    if( t == 0 || D == 0 )
    {
        return 0.0;
    }

    if( getr0() > getrsink() )
        return flux_abs_Ll( t );
    else
        return flux_abs_Lr( t );
}

/* Function returns flux at absorbing bounday a. */
Real GreensFunction1DAbsSinkAbs::flux_leavea(Real t) const
{
    if( t == 0 || D == 0 )
    {
        return 0.0;
    }

    if( getr0() < getrsink() )
        return flux_abs_Ll( t );
    else
        return flux_abs_Lr( t );
}

/* Calculates the total probability flux leaving the domain at time t
   This is simply the negative of the time derivative of the survival prob.
   at time t [-dS(t')/dt' for t'=t]. */
Real GreensFunction1DAbsSinkAbs::flux_tot(Real t) const
{
    const Real D( getD() );
    const Real k( getk() ); 
    const Real Lr( getLr() );
    const Real Ll( getLl() );
    const Real L0( getL0() );

    const Real L( Lr + Ll );
    const Real LlpL0( Ll + L0 );
    const Real LrmL0( Lr - L0 );

    Real sum = 0, term_n = 0, prev_term = 0;
    Real root_n, root_n2, term1, term2, numerator;
    int n = 1;
    do
    {
        if ( n >= MAX_TERMS )
        {
            std::cerr << " Too many terms needed for GF1DSink::flux_tot. N: "
                      << n << std::endl;
            break;
        }
        prev_term = term_n;
    
        root_n = get_root(n);
        root_n2 = gsl_pow_2( root_n );

        term1 = sin( root_n * L ) - sin( root_n * LrmL0 ) - sin( root_n * LlpL0 );
        term2 = sin( root_n * Lr ) - sin( root_n * L0 ) - sin( root_n * LrmL0 );
        
        numerator = D * term1 + k * sin( root_n * Ll ) * term2 / root_n;

        term_n = root_n2 * Greens_fn(t, root_n, root_n2) * numerator;
        sum += term_n;

        n++;
    }
    while ( fabs(term_n/sum) > EPSILON  ||
	        fabs(prev_term/sum) > EPSILON ||
	        n <= MIN_TERMS );

    return  2.0 * D * sum;
}

/* Flux leaving throught absorbing boundary from sub-domain containing r0. */
Real GreensFunction1DAbsSinkAbs::flux_abs_Lr(Real t) const
{
    const Real D( getD() );
    const Real k( getk() );    
    const Real Ll( getLl() );
    const Real L0( getL0() );
    const Real LlpL0( Ll + L0 );
    Real root_n;
    
    Real sum = 0, numerator, term_n = 0, prev_term = 0;
    int n = 1;
    do
    {
        if ( n >= MAX_TERMS )
        {
            std::cerr << " Too many terms needed for GF1DSink::flux_abs_Lr. N: "
                      << n << std::endl;
            break;
        }
        prev_term = term_n;
    
        root_n = get_root(n);
        
        numerator = k * sin( root_n * Ll ) * sin ( root_n * L0 ) + D * root_n * sin( root_n * LlpL0 );
        numerator *= root_n;

        term_n = Greens_fn(t, root_n, gsl_pow_2( root_n ) ) * numerator;         
        sum += term_n;

        n++;
    }
    while ( fabs(term_n/sum) > EPSILON  ||
            fabs(prev_term/sum) > EPSILON ||
            n <= MIN_TERMS );

    return - D * 2 * sum;
}

/* Flux leaving throught absorbing boundary from sub-domain not containing r0. */
Real GreensFunction1DAbsSinkAbs::flux_abs_Ll(Real t) const
{
    const Real D2( gsl_pow_2( getD() ) );
    const Real LrmL0( getLr() - getL0() );

    Real root_n, root_n2;
    
    Real sum = 0, numerator = 0, term_n = 0, prev_term = 0;
    int n = 1;
    do
    {
        if ( n >= MAX_TERMS )
        {
            std::cerr << " Too many terms needed for GF1DSink::flux_abs_Ll. N: "
                      << n << std::endl;
            break;
        }
        prev_term = term_n;
    
        root_n = get_root(n);
        root_n2 = gsl_pow_2( root_n );        
        
        numerator = root_n2 * sin( root_n * LrmL0 );

        term_n = Greens_fn(t, root_n, root_n2 ) * numerator;
        sum += term_n;

        n++;
    }
    while ( fabs(term_n/sum) > EPSILON  ||
            fabs(prev_term/sum) > EPSILON ||
            n <= MIN_TERMS );
    
    return 2 * D2 * sum;
}

/* Calculates the probability flux leaving the domain through the sink
   at time t */
Real GreensFunction1DAbsSinkAbs::flux_sink(Real t) const
{
    return getk() * prob_r(getrsink(), t);
}

/* Determine which event has occured, an escape or a reaction. Based on the
   fluxes through the boundaries and the sink at the given time. */
GreensFunction1DAbsSinkAbs::EventKind 
GreensFunction1DAbsSinkAbs::drawEventType( Real rnd, Real t ) const
{
    THROW_UNLESS( std::invalid_argument, rnd < 1.0 && rnd >= 0.0 );
    THROW_UNLESS( std::invalid_argument, t > 0.0 );

    const Real a( geta() );
    const Real sigma( getsigma() );
    const Real r0( getr0() );
    const Real L( a - sigma );

    // if the radiative boundary is impermeable (k==0) or
    // the particle is at one the absorbing boundary (sigma or a) => IV_ESCAPE event
    if ( k == 0 || fabs(a - r0) < EPSILON * L || fabs(sigma - r0) < EPSILON * L)
    {
	    return IV_ESCAPE;
    }

    /* The event is sampled from the flux ratios.
       Two possiblities:
       (1) Leave through left or right boundary - IV_ESCAPE
       (2) Leave through sink - IV_REACTION
    */
    rnd *= flux_tot( t );
    const Real p_sink( flux_sink( t ) );
    if (rnd < p_sink)
    {
      return IV_REACTION;
    }
    else
    {
      return IV_ESCAPE;
    }
}

/* This function is needed to cast the math. form of the function
   into the form needed by the GSL root solver. */
Real GreensFunction1DAbsSinkAbs::drawT_f (Real t, void *p)
{
    // casts p to type 'struct drawT_params *'
    struct drawT_params *params = (struct drawT_params *)p;
    Real Xn = 0, exponent = 0;
    Real prefactor = params->prefactor;
    int terms = params->terms;

    Real sum = 0, term = 0, prev_term = 0;
    int n = 0;
    do
    {
        if ( n >= terms )
        {
            std::cerr << "Too many terms needed for GF1DSink::DrawTime. N: "
                      << n << ", conv arg: " << exponent * t <<  std::endl;
            break;
        }
        prev_term = term;
        
        Xn = params->Xn[n];
        exponent = params->exponent[n];
        term = Xn * exp(exponent * t);
	
        sum += term;
        n++;
    }
    while (fabs(term/sum) > EPSILON*1.0 ||
           fabs(prev_term/sum) > EPSILON*1.0 ||
           n <= MIN_TERMS );

    // find the intersection with the random number
    return params->rnd - prefactor*sum;
}

// Draws the first passage time from the survival probability,
// using an assistance function drawT_f that casts the math. function
// into the form needed by the GSL root solver.
Real GreensFunction1DAbsSinkAbs::drawTime(Real rnd) const
{
    THROW_UNLESS( std::invalid_argument, 0.0 <= rnd && rnd < 1.0 );
  
    const Real a( geta() );    
    const Real r0( getr0() );
    const Real D( getD() );
    const Real k( getk() ); 
    const Real Lr( getLr() );
    const Real Ll( getLl() );
    const Real L0( getL0() );
    const Real L( getLr() + getLl() );

    Real root_n, Xn, exponent, term1, term2, numerator;

    if ( D == 0.0 || L == INFINITY )
    {
	    return INFINITY;
    }

    if ( rnd <= EPSILON || L < 0.0 || fabs( a - r0 ) < EPSILON * L )
    {
	    return 0.0;
    }

    // the structure to store the numbers to calculate the numbers for 1-S
    struct drawT_params parameters;

    /* produce the coefficients and the terms in the exponent and put them
       in the params structure. This is not very efficient at this point,
       coefficients should be calculated on demand->TODO or not TODO? */
    for (int n = 0; n < MAX_TERMS; n++)
    {
        root_n = get_root( n + 1 );	

        term1 = sin( root_n * L ) - sin( root_n * (Lr - L0) ) - sin( root_n * (Ll + L0) );
        term2 = sin( root_n * Lr ) - sin( root_n * L0 ) - sin( root_n * (Lr - L0) );
        
        numerator = D * term1 + k * sin( root_n * Ll ) * term2 / root_n;
        Xn = numerator / p_denominator( root_n );
		      
        exponent = - D * gsl_pow_2( root_n );

        // store the coefficients in the structure
        parameters.Xn[n] = Xn;
	    // also store the values for the exponent
        parameters.exponent[n] = exponent;
    }
    
    // the prefactor of the sum is also different in case of drift<>0 :
    parameters.prefactor = 2.0;
    
    // store the random number for the probability
    parameters.rnd = rnd;
    // store the number of terms used
    parameters.terms = MAX_TERMS;
    parameters.tscale = t_scale;

    // Define the function for the rootfinder
    gsl_function F;
    F.function = &GreensFunction1DAbsSinkAbs::drawT_f;
    F.params = &parameters;

    // Find a good interval to determine the first passage time in
    // get the distance to absorbing boundary (disregard Sink)
    Real t_guess;
    const Real dist( std::min(Lr, Ll) );
    
    t_guess = dist * dist / ( 2.0 * D );
    t_guess *= .1;

    Real value( GSL_FN_EVAL( &F, t_guess ) );
    Real low( t_guess );
    Real high( t_guess );

    // scale the interval around the guess such that the function straddles
    if( value < 0.0 )
    {
        // if the guess was too low
        do
        {
            // keep increasing the upper boundary until the
            // function straddles
            high *= 10;
            value = GSL_FN_EVAL( &F, high );
            
            if( fabs( high ) >= t_guess * 1e6 )
            {
                std::cerr << "GF1DSink::drawTime Couldn't adjust high. F("
                          << high << ") = " << value << std::endl;
                throw std::exception();
            }
        }
        while ( value <= 0.0 );
    }
    else
    {
        // if the guess was too high
	    // initialize with 2 so the test below survives the first
	    // iteration
	    Real value_prev( 2 );
	    do
	    {
	        if( fabs( low ) <= t_guess * 1e-6 ||
	            fabs(value-value_prev) < EPSILON*1.0 )
	        {
		    std::cerr << "GF1DSink::drawTime Couldn't adjust low. F(" << low << ") = "
		              << value << " t_guess: " << t_guess << " diff: "
		              << (value - value_prev) << " value: " << value
		              << " value_prev: " << value_prev << " rnd: "
		              << rnd << std::endl;
		    return low;
	        }
	        value_prev = value;
	        // keep decreasing the lower boundary until the function straddles
	        low *= 0.1;
	        // get the accompanying value
	        value = GSL_FN_EVAL( &F, low );
	    }
	    while ( value >= 0.0 );
    }

    // find the intersection on the y-axis between the random number and
    // the function
    
    // define a new solver type brent
    const gsl_root_fsolver_type* solverType( gsl_root_fsolver_brent );
    // make a new solver instance
    gsl_root_fsolver* solver( gsl_root_fsolver_alloc( solverType ) );
    const Real t( findRoot( F, solver, low, high, t_scale*EPSILON, EPSILON,
                            "GreensFunction1DAnsSinkAbs::drawTime" ) );
    // return the drawn time
    return t;
}

Real GreensFunction1DAbsSinkAbs::drawR_f(Real rr, void *p)
{
    // casts p to type 'struct drawR_params *'
    struct drawR_params *params = (struct drawR_params *)p;
    int  max_terms = params->terms;
    
    Real sum = 0, term_n = 0, prev_term = 0;
    
    /* Determine in which part of the domain rr lies, and
       thus which domain function to use. */
    static Real (*p_int_r)
        (Real const&, Real const&, drawR_params const&) = NULL;

    if( rr <= 0 )
        p_int_r = &GreensFunction1DAbsSinkAbs::p_int_r_leftdomain;
    else if( rr < params->L0 )
            p_int_r = &GreensFunction1DAbsSinkAbs::p_int_r_rightdomainA;
        else
            p_int_r = &GreensFunction1DAbsSinkAbs::p_int_r_rightdomainB;
    
    int n = 0;
    do
    {
	    if ( n >= max_terms )
	    {
	        std::cerr << "GF1DSink: Too many terms needed for DrawR. N: "
	                  << n << std::endl;
	        break;
	    }
	    prev_term = term_n;

	    term_n = params->exp_and_denominator[n] * (*p_int_r)
            ( rr, params->root_n[n], *params );

	    sum += term_n;
	    n++;
    }
    while (fabs(term_n/sum) > EPSILON*1.0 ||
	       fabs(prev_term/sum) > EPSILON*1.0 ||
           n <= MIN_TERMS );
    
    // Find the intersection with the random number
    return params->prefactor * sum - params->rnd;
}

//Integrated Greens function for rr part of [-Ll, 0]
Real GreensFunction1DAbsSinkAbs::p_int_r_leftdomain(Real const& rr, 
                                                    Real const& root_n,
                                                    drawR_params const& params)
{
    const Real LrmL0( params.Lr - params.L0 );
    const Real Llprr( params.Ll + rr );
   
    return params.D * sin( root_n * LrmL0 ) * ( cos( root_n * Llprr ) - 1.0 );
}

//Integrated Greens function for rr part of (0, L0]
Real GreensFunction1DAbsSinkAbs::p_int_r_rightdomainA(Real const& rr, 
                                                      Real const& root_n, 
                                                      drawR_params const& params)
{
    const Real LrmL0( params.Lr - params.L0 );
    const Real Llprr( params.Ll + rr );
    const Real root_n_rr( root_n * rr );
    
    const Real temp( params.D * ( cos( root_n * Llprr ) - 1.0 ) + 
            params.k / root_n * ( cos( root_n_rr ) - 1 ) * sin( root_n * params.Ll ) );
    
    return sin( root_n * LrmL0 ) * temp;
}

//Integrated Greens function for rr part of (L0, Lr]
Real GreensFunction1DAbsSinkAbs::p_int_r_rightdomainB(Real const& rr, 
                                                      Real const& root_n, 
                                                      drawR_params const& params)
{
    const Real Lr( params.Lr );
    const Real Ll( params.Ll );
    const Real L0( params.L0 );
    const Real LrmL0( Lr - L0 );
    const Real Lrmrr( Lr - rr );
    const Real L( Lr + Ll );
    const Real LlpL0( Ll + L0 );
                   
    const Real term1( sin( root_n * L ) - sin( root_n * LrmL0 ) - 
            sin( root_n * LlpL0 ) * cos( root_n * Lrmrr ) );
            
    const Real term2( sin( root_n * Lr ) - sin( root_n * LrmL0 ) -
            sin( root_n * L0 ) * cos( root_n * Lrmrr ));
        
    return params.D * term1 + params.k * sin( root_n * Ll ) * term2 / root_n;
}

Real GreensFunction1DAbsSinkAbs::drawR(Real rnd, Real t) const
{
    THROW_UNLESS( std::invalid_argument, 0.0 <= rnd && rnd <= 1.0 );
    THROW_UNLESS( std::invalid_argument, t >= 0.0 );
    
    const Real D( getD() );
    const Real Lr( getLr() );
    const Real Ll( getLl() );
    const Real rsink( getrsink() );
    const Real L( Lr + Ll );

    if (t == 0.0 || (D == 0.0 && v == 0.0) )
    {
	    // the trivial case
	    return r0;
    }

    if ( L < 0.0 )
    {
	    // if the domain had zero size
	    return 0.0;
    }

    if ( rnd <= EPSILON )
    {
        return -Ll;        
    }

    if( rnd >= ( 1 - EPSILON ) )
    {
        return Lr;
    }

    // the structure to store the numbers to calculate r.
    struct drawR_params parameters;
    const Real S( p_survival(t) );
    Real root_n, root_n2;

    assert(S >= 0.0);

    /* produce the coefficients and the terms in the exponent and put them
       in the params structure //TODO: MAX_TERM -> CONVERGENCE_TERMS */
    for (int n = 0; n < MAX_TERMS; n++)
    {
	    root_n = get_root(n+1);
       	root_n2 = gsl_pow_2( root_n );

        // also store the exponent and denominator.
        parameters.exp_and_denominator[n] = Greens_fn( t, root_n, root_n2 );
    }

    //Store pointer to the rootList in params.
    parameters.rootListPtr = get_rootListPtr();
    
    // The normalization.
    parameters.prefactor = 2.0 / S;

    // and the rest of the parameters
    parameters.rnd = rnd;
    parameters.L0 = getL0();
    parameters.Lr = getLr();
    parameters.Ll = getLl();
    parameters.D = getD();
    parameters.k = getk();

    // store the number of terms used
    parameters.terms = MAX_TERMS;

    // define gsl function for rootfinder
    gsl_function F;
    F.function = &drawR_f;
    F.params = &parameters;

    // define a new solver type brent
    const gsl_root_fsolver_type* solverType( gsl_root_fsolver_brent );
    // make a new solver instance
    // TODO: incl typecast?
    gsl_root_fsolver* solver( gsl_root_fsolver_alloc( solverType ) );
    Real rr( findRoot( F, solver, -Ll, Lr, EPSILON*L, EPSILON,
                            "GreensFunction1DRadAbs::drawR" ) );

    // Convert the position rr to 'world' coordinates and return it.
    if( getr0() < rsink )
        return rsink - rr;
    else
        return rsink + rr;
}

std::string GreensFunction1DAbsSinkAbs::dump() const
{
    std::ostringstream ss;
    ss << "D = " << getD() << ", sigma = " << getsigma() <<
        ", a = " << geta() <<
        ", r0 = " << getr0() <<
        ", rsink = " << getrsink() <<
        ", k = " << getk() << std::endl;
    return ss.str();
}

