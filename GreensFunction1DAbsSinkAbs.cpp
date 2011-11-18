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
#include "GreensFunction1DAbsSinkAbs.hpp"

// This is the appropriate definition of the function defining
// the roots of our Green's functions in GSL.
// Later needed by the rootfinder.
Real GreensFunction1DAbsSinkAbs::root_f (Real x, void *p)
{
    // casts the void to the struct pointer
    struct root_f_params *params = (struct root_f_params *)p;
    const Real Lm_L = (params->Lm_L);
    const Real h = (params->h);

    // L   = Lr + Ll
    // h    = L * k / (2 * D)
    // L_Lm = Lr + Ll / Lr - Ll
    // x    = q * L

    return x * sin(x) + h * ( cos(x * Lm_L) - cos(x) );
}

// Calculates the roots of tan(x*a)=-x/h
Real GreensFunction1DAbsSinkAbs::root_n(int const& n) const
{
    const Real L( getLr() + getLl() );
    const Real Lm( getLr() - getLl() );
    const Real h( getk() * L / ( 2 * getD() ) );
    
    Real upper, lower;  

    //h sets how strong second term in root function is pronounced. 
    if ( h < 1 )
    {
	    // 1E-10 to make sure that he doesn't include the transition 
	    lower = (n-1)*M_PI + 1E-10;
	    // (asymptotic) from infinity to -infinity
	    upper =  n   *M_PI - 1E-10;
    }
    else
    {
        //TODO: more complicated scheme needed.
        THROW( h > 1 )
    }

    gsl_function F;
    struct root_f_params params = { Lm / L, h };
     
    F.function = &GreensFunction1DAbsSinkAbs::root_f;
    F.params = &params;

    // define a new solver type brent
    const gsl_root_fsolver_type* solverType( gsl_root_fsolver_brent );

    // make a new solver instance
    // TODO: incl typecast?
    gsl_root_fsolver* solver( gsl_root_fsolver_alloc( solverType ) );
    // get the root = run the rootsolver
    const Real root( findRoot( F, solver, lower, upper, 1.0*EPSILON, EPSILON,
                            "GreensFunction1DAbsSinkAbs::root_f" ) );
    gsl_root_fsolver_free( solver );
    
    return root/L;
}

/* Denominator of the greensfunction. */
Real GreensFunction1DAbsSinkAbs::p_denominator(Real const& root_n) const
{
    const Real D( getD() );
    const Real k( getk() );    
    const Real Lm( getLr() - getLl() );
    const Real L( getLr() + getLl() );
    
    const Real term1( root_n * L * cos( root_n * L ) + sin( root_n * L );
    const Real term2( L * sin( root_n * L ) + Lm * sin( root_n * L );
    
    return D * term1 + k / 2. * term2;
}


// Calculates the probability of finding the particle inside the domain
// at time t, the survival probability.
Real GreensFunction1DAbsSinkAbs::p_survival(Real const& t) const
{
    THROW_UNLESS( std::invalid_argument, t >= 0.0 );
  
    const Real D ( getD() );
    const Real k ( getk() ); 
    const Real L ( getLr() + getLl() );
    const Real Lr( getLr() );
    const Real Ll( getLl() );
    const Real L0( getL0() );

    if (t == 0.0 || (D == 0.0 && v == 0.0) )
    {
	    //particle can't escape.
	    return 1.0;
    }

    Real root_n;
    Real sum = 0, num = 0, num_prev = 0;
    int n = 1;
    do
    {
        if ( n >= MAX_TERMS )
	    {
	        std::cerr << " Too many terms needed for GF1DSink::p_survival. N: "
	              << n << std::endl;
	        break;
	    }
	    num_prev = num;
    
        root_n = root_n(n);
        const Real term1( sin( root_n * L ) - sin( root_n * (Lr - L0) ) - sin( root_n * (Ll + L0) ) );
        const Real term2( sin( root_n * Lr ) - sin( root_n * L0 ) - sin( root_n * (Lr - L0) ) );
        
        num = ( D * term1 + k * sin( root_n * Ll ) * term2 / root_n );
	    sum += exp( - D * gsl_pow_2( root_n ) * t ) * num / p_denominator( root_n );
	    n++;
    }
    while ( fabs(term/sum) > EPSILON  ||
	        fabs(term_prev/sum) > EPSILON ||
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
    const Real L( getLr() + getLl() );
    const Real Lr( getLr() );
    const Real Ll( getLl() );
    Real L0( getL0() );
    Real rr( r - getrsink() ) ;

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
           
    /* If r is on the right side of the sink, use right the solution for the right domain.
       Otherwise, use the solution for the left domain. */
    if( rr >= 0 )
    {
        /* Check if r is left or right of the starting position r0.
           If so, interchange r with r0. */
        if( rr < L0 )
        {
            const Real temp( L0 );
            rr = L0;
            L0 = temp;
        }    

        do
        {   
	        if ( n >= MAX_TERMS )
	        {
	            std::cerr << " Too many terms needed for GF1DSink::prob_r. N: "
	              << n << std::endl;
	            break;
	        }
	        term_prev = term_n;

	        root_n = root_n(n);

            numerator =  D * root_n * sin( root_n * (Ll + L0) ) + k * sin( root_n * Ll ) * sin( root_n * x0);
            numerator *= sin( root_n * (Lr - rr) );       
                
            term_n = exp( - D * gsl_pow_2( root_n ) * t ) * numerator / p_denominator( root_n );
                
    	    sum += term_n;

    	    n++;
       }
       while ( fabs(term/sum) > EPSILON*PDENS_TYPICAL || 
	           fabs(prev_term/sum) > EPSILON*PDENS_TYPICAL ||
	           n <= MIN_TERMS );
	}
    else
    {
        do
        {   
	        if ( n >= MAX_TERMS )
	        {
	            std::cerr << " Too many terms needed for GF1DSink::prob_r. N: "
	              << n << std::endl;
	            break;
	        }
	        term_prev = term_n;

	        root_n = root_n(n);

            numerator =  D * root_n * sin( root_n * (Ll + rr) ) * sin( root_n * (Lr - L0) ); 
            
            term_n = exp( - D * gsl_pow_2( root_n ) * t ) * numerator / p_denominator( root_n );
               
    	    sum += term_n;

    	    n++;
       }
       while ( fabs(term/sum) > EPSILON*PDENS_TYPICAL || 
	           fabs(prev_term/sum) > EPSILON*PDENS_TYPICAL ||
	           n <= MIN_TERMS );
    }
    
    return - 2.0 * sum;
}

// Calculates the probability density of finding the particle at location z at
// timepoint t, given that the particle is still in the domain.
Real GreensFunction1DAbsSinkAbs::calcpcum(Real r, Real t) const
{
    return prob_r(r, t)/p_survival(t);
}

// Calculates the total probability flux leaving the domain at time t
// This is simply the negative of the time derivative of the survival prob.
// at time t [-dS(t')/dt' for t'=t].
Real GreensFunction1DAbsSinkAbs::flux_tot(Real t) const
{
    Real root_n;
    const Real D( getD() );
    const Real k( getk() ); 
    const Real L( getLr() + getLl() );
    const Real Lr( getLr() );
    const Real Ll( getLl() );
    const Real L0( getL0() );
    const Real rr( r - getrsink() ) ;

    double sum = 0, term = 0, prev_term = 0;
    int n = 1;
    do
    {
        if ( n >= MAX_TERMS )
	    {
	        std::cerr << " Too many terms needed for GF1DSink::flux_tot. N: "
	              << n << std::endl;
	        break;
	    }
	    num_prev = num;
    
        root_n = root_n(n);
        const Real term1( sin( root_n * L ) - sin( root_n * (Lr - L0) ) - sin( root_n * (Ll + L0) ) );
        const Real term2( sin( root_n * Lr ) - sin( root_n * L0 ) - sin( root_n * (Lr - L0) ) );
        
        num = ( D * term1 + k * sin( root_n * Ll ) * term2 / root_n );
	    sum += gsl_pow_2( root_n ) * exp( - D * gsl_pow_2( root_n ) * t ) * num / p_denominator( root_n );
	    n++;
    }
    while ( fabs(term/sum) > EPSILON  ||
	        fabs(term_prev/sum) > EPSILON ||
	        n <= MIN_TERMS );

    return  2.0 * D * sum;
}

/* Calculates the probability flux leaving the domain through the absorbing
   boundaries at the left and right side at time t. */
real_pair GreensFunction1DAbsSinkAbs::flux_abs(Real t) const
{
    const Real D( getD() );
    const Real k( getk() );    
    const Real Lr_L0( getLr() - getL0() );
    const Real LlpL0( getLr() + getL0() );
    const Real Ll( getLl() );
    const Real L0( getL0() );
    Real root_n;
    
    Real sum = 0, term = 0, prev_term = 0, prefac = 0;
    int n = 1;
    do
    {
        if ( n >= MAX_TERMS )
	    {
	        std::cerr << " Too many terms needed for GF1DSink::flux_abs. N: "
	              << n << std::endl;
	        break;
	    }
	    num_prev = num_left;
    
        root_n = root_n(n);
        root_n2 = gsl_pow_2( root_n );        
        
        num_left = root_n2 * sin( root_n * ( Lr_L0 ) );

        num_right = k * sin( root_n * Ll ) * sin ( root_n * L0 ) + D * root_n * sin( root_n * LlpL0 );
        
        prefac = exp( - D * root_n2 * t ) / p_denominator( root_n );
        
	    sum_left  += prefac * num_left;
	    sum_right += prefac * num_right; 
	    n++;
    }
    while ( fabs(term/sum_left) > EPSILON  ||
	        fabs(term_prev/sum_left) > EPSILON ||
	        n <= MIN_TERMS );

    return real_pair( 2 * D * D * sum_left,  - D * 2 * root_n * sum_right);
}

// Calculates the probability flux leaving the domain through the sink
// at time t
Real GreensFunction1DAbsSinkAbs::flux_sink(Real t) const
{
    return getk() * prob_r(getrsink(), t);
}

/* Determine which event has occured, an escape or a reaction. Based on the
   fluxes through the boundaries and the sink at the given time. Beware: if t is not a
   first passage time you still get an answer! */
GreensFunction1DAbsSinkAbs::EventKind
GreensFunction1DAbsSinkAbs::drawEventType( Real rnd, Real t )
const
{
    THROW_UNLESS( std::invalid_argument, rnd < 1.0 && rnd >= 0.0 );
    THROW_UNLESS( std::invalid_argument, t > 0.0 );
    // if t=0 nothing has happened => no event

    const Real  a( geta() );
    const Real  L( geta() - getsigma() );
    const Real r0( getr0() );

    // if the radiative boundary is impermeable (k==0) or
    // the particle is at one the absorbing boundary (sigma or a) => IV_ESCAPE event
    if ( k == 0 || fabs(a - r0) < EPSILON * L || fabs(sigma - r0) < EPSILON * L)
    {
	    return IV_ESCAPE;
    }

    /* The event is sampled from the flux ratios.
       Three possiblities:
       (1) Leave through left  boundary - IV_ESCAPE_L
       (2) Leave through right boundary - IV_ESCAPE_R
       (3) Leave through sink - IV_REACTION
    */
    rnd *= flux_tot( t );
    //Get fluxes through left and right boundaries.
    const real_pair flux_abs_left_right( flux_abs( t ) );
    
    Real p_accumulate( flux_abs_left_right.first );
    if (rnd < p_accumulate)
    {
	    return IV_ESCAPE_L;
    }
    else
    {
        p_accumulate += flux_abs_left_right.second;
        if( rnd < p_accumulate )
        {
	        return IV_ESCAPE_R;
	    }
	    else
	    {
	        return IV_REACTION;
	    }
    }
}

// This function is needed to cast the math. form of the function
// into the form needed by the GSL root solver.
Real GreensFunction1DAbsSinkAbs::drawT_f (Real t, void *p)
{
    // casts p to type 'struct drawT_params *'
    struct drawT_params *params = (struct drawT_params *)p;
    Real Xn, exponent;
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
    return 1.0 - prefactor*sum - params->rnd;
}

// Draws the first passage time from the survival probability,
// using an assistance function drawT_f that casts the math. function
// into the form needed by the GSL root solver.
Real GreensFunction1DAbsSinkAbs::drawTime(Real rnd) const
{
    THROW_UNLESS( std::invalid_argument, 0.0 <= rnd && rnd < 1.0 );
  
    const Real sigma( getsigma() );
    const Real a( geta() );    
    const Real r0( getr0() );
    const Real D( getD() );
    const Real k( getk() ); 
    const Real L( getLr() + getLl() );
    const Real Lr( getLr() );
    const Real Ll( getLl() );
    const Real L0( getL0() );
    const Real rr( r - getrsink() ) ;

    Real root_n;    

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
    // some temporary variables
    double root_n = 0;
    double root_n2, root_n_r0_s, root_n_L, h_root_n;
    double Xn, exponent, prefactor;

    // produce the coefficients and the terms in the exponent and put them
    // in the params structure. This is not very efficient at this point,
    // coefficients should be calculated on demand->TODO or not TODO?
    for (int n=0; n<MAX_TERMS; n++)
    {
	    root_n = root_n( n + 1 );	
		
        const Real term1( sin( root_n * L ) - sin( root_n * (Lr - L0) ) - sin( root_n * (Ll + L0) ) );
        const Real term2( sin( root_n * Lr ) - sin( root_n * L0 ) - sin( root_n * (Lr - L0) ) );
        
        const Real numerator( D * term1 + k * sin( root_n * Ll ) * term2 / root_n );
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
.    
    Real sum = 0, term = 0, prev_term = 0;
    Real root_n, exp_and_denominator;
    
    /* Determine in which part of the domain rr lies, and
       thus which domain function to use. */
    Real (*numerator_int_f)(Real const&, Real const&, Real const&) = NULL;
    if( rr <= 0 )
        numerator_int_f = &GreensFunction1DAbsSinkAbs::num_int_r_leftdomain;
    else if( rr < L0 )
            numerator_int_f = &GreensFunction1DAbsSinkAbs::num_int_r_rightdomainA;
        else
            numerator_int_f = &GreensFunction1DAbsSinkAbs::num_int_r_rightdomainB;
    
    int n = 0;
    do
    {
	    if ( n >= max_terms )
	    {
	        std::cerr << "GF1DSink: Too many terms needed for DrawR. N: "
	                  << n << std::endl;
	        break;
	    }
	    prev_term = term;

	    exp_and_denominator = params->exp_and_denominator[n];	    
	    term = exp_and_denominator * (*numerator_int_f)( rr, t, root_n[n] );

	    sum += term;
	    n++;
    }
    while (fabs(term/sum) > EPSILON*1.0 ||
	       fabs(prev_term/sum) > EPSILON*1.0 ||
           n <= MIN_TERMS );

    // Find the intersection with the random number
    return params.prefactor * sum - params->rnd;
}

//Integrated Greens function for rr part of [-Ll, 0]
inline Real GreensFunction1DAbsSinkAbs::num_int_r_leftdomain(Real const& rr, Real const& t, Real const& root_n) const
{
    const Real Ll_L0( getLl() - getL0() );
    const Real Ll_rr( getLl() + rr );
    
    return getD() * sin( root_n * Ll_L0 ) * ( cos( root_n * Ll_rr ) - 1 )
}

//Integrated Greens function for rr part of (0, L0]
inline Real GreensFunction1DAbsSinkAbs::num_int_r_rightdomainA(Real const& rr, Real const& t, Real const& root_n) const
{
    const Real Lr_L0( getLr() - getL0() );
    const Real Ll_rr( getLl() + rr );
    const Real root_n_rr( root_n * rr );
    
    const Real temp( getD() * ( cos( root_n * Ll_rr ) - 1 ) + 
            getk() / root_n * ( cos( root_n_rr ) - 1 ) * sin( root_n * Ll ) );
    
    return sin( root_n * Lr_L0 ) * temp;
}

//Integrated Greens function for rr part of (L0, Lr]
inline Real GreensFunction1DAbsSinkAbs::num_int_r_rightdomainB(Real const& rr, Real const& t, Real const& root_n) const
{
    const Real L( getLr() + getLl() );
    const Real Lr( getLr() );
    const Real Ll( getLl() );
    const Real L0( getL0() );
    const Real Lr_L0( Lr - L0 );
    const Real Lr_rr( Lr - rr );
                   
    const Real term1( sin( root_n * L ) - sin( root_n * Lr_L0 ) - 
            sin( root_n * (Ll + L0) ) * cos( root_n * Lr_rr ) );
            
    const Real term2( sin( root_n * Lr ) - sin( root_n * Lr_L0 ) -
            sin( root_n * L0 ) * cos( root_n * Lr_rr ));
        
    return getD() * term1 + getk() * sin( root_n * Ll ) * term2 / root_n;
}

Real GreensFunction1DAbsSinkAbs::drawR(Real rnd, Real t) const
{
    THROW_UNLESS( std::invalid_argument, 0.0 <= rnd && rnd < 1.0 );
    THROW_UNLESS( std::invalid_argument, t >= 0.0 );
    
    const Real D( getD() );
    const Real L( getLr() + getLl() );
    const Real Lr( getLr() );
    const Real Ll( getLl() );
    const Real L0( getL0() );

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

    // the structure to store the numbers to calculate r.
    struct drawR_params parameters;
    const Real S( p_survival(t) );
    Real root_n = 0;

    assert(S >= 0.0);

    // produce the coefficients and the terms in the exponent and put them
    // in the params structure //TODO: MAX_TERM -> CONVERGENCE_TERMS
    for (int n = 0; n < MAX_TERMS; n++)
    {
	    root_n = root_n(n+1);
       	root_n2 = root_n * root_n;
      	// store the roots in the structure
	    parameters.root_n[n] = root_n;
        // also store the exponent and denominator.
        parameters.exp_and_denominator[n] = exp( - D * root_n2 * t ) / p_denominator( root_n );
        // and the normalization.
        parameters.prefactor = 2.0 / S;
    }
    
    // define gsl function for rootfinder
    gsl_function F;
    F.function = &GreensFunction1DRadAbs::drawR_f;
    F.params = &parameters;

    // store the random number for the probability
    parameters.rnd = rnd;
    // store the number of terms used
    parameters.terms = MAX_TERMS;

    // define a new solver type brent
    const gsl_root_fsolver_type* solverType( gsl_root_fsolver_brent );
    // make a new solver instance
    // TODO: incl typecast?
    gsl_root_fsolver* solver( gsl_root_fsolver_alloc( solverType ) );
    Real r( findRoot( F, solver, -Ll, Lr, EPSILON*L, EPSILON,
                            "GreensFunction1DRadAbs::drawR" ) );

    // Convert the position to 'world' coordinates and return it.
    
    return r;
}

std::string GreensFunction1DAbsSinkAbs::dump() const
{
    std::ostringstream ss;
    ss << "D = " << this->getD() << ", sigma = " << this->getsigma() <<
        ", a = " << this->geta() <<
        ", r0 = " << this->getr0() <<
        ", rsink = " << this->getrsink() <<
        ", k = " << this->getk() << std::endl;
    return ss.str();
}

