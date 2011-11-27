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
#include "funcSum.hpp"
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
void GreensFunction1DAbsSinkAbs::calculate_n_roots(uint const& n)
{
    const Real Lr( getLr() );
    const Real Ll( getLl() );
    const Real L( Lr + Lr );
    const Real Lm( Lr - Ll );
    const Real Lm_L( Lm / L );
    const Real h( getk() * L / ( 2 * getD() ) );
    
    Real root_i;
    uint i = 0;

    real_pair lower_upper_pair;

    /* set params needed for get lower and upper function */
    lo_up_params.h = h;
    lo_up_params.Lm_L = Lm_L;
    lo_up_params.long_half_wave = std::max( L/Lr * M_PI, L/Ll * M_PI );
    lo_up_params.short_half_wave = std::min( L/Lr * M_PI, L/Ll * M_PI );

    /* Define the root function. */
    gsl_function F;
    struct root_f_params params = { Lm_L, h };
     
    F.function = &GreensFunction1DAbsSinkAbs::root_f;
    F.params = &params;

    /* define and ad a new solver type brent */
    const gsl_root_fsolver_type* solverType( gsl_root_fsolver_brent );
    gsl_root_fsolver* solver( gsl_root_fsolver_alloc( solverType ) );

    /* Find all the roots up to the nth*/
    i = rootList_size();

    /* If rootList is empty, add first root: 0. Will be removed later. */
    uint n_adapt( n );
    if(i == 0)
    {
        /* Stores pointers to last root with long/short wavelength. */
        lo_up_params.last_long_root = 0;
        lo_up_params.last_short_root = 0;
        ad_to_rootList( 0.0 );
        n_adapt++;
    }

    while(i++ < n_adapt)
    {
        lower_upper_pair = get_lower_and_upper( lo_up_params );

        root_i = findRoot( F, solver, lower_upper_pair.first, lower_upper_pair.second, 
                           1.0*EPSILON, EPSILON, "GreensFunction1DAbsSinkAbs::root_f" );

        assert( root_i > get_last_root() - EPSILON );

        ad_to_rootList( root_i );
    }

    /* Remove first root */
    if( get_root( 0 ) == 0.0 )
    {
        remove_from_rootList( 0 );
        /* correct pointers for removal */
        if( lo_up_params.last_long_root > 0 ) 
            lo_up_params.last_long_root--;
        if( lo_up_params.last_short_root > 0 ) 
            lo_up_params.last_short_root--;
    }

    gsl_root_fsolver_free( solver );
}


/* returns two points on the x-axis which straddle the next root. */
std::pair<Real, Real> GreensFunction1DAbsSinkAbs::get_lower_and_upper( lower_upper_params& params )
{
    const Real root_n( get_last_root() );
    const Real safety( .75 );

    Real lower, upper, next_root_est, left_offset, right_offset;

    if( params.h / root_n < 1 )
    {
        right_offset = M_PI;
        next_root_est = root_n + M_PI;
    }
    else
    {
        const Real next_root_long( get_root( params.last_long_root ) + params.long_half_wave );
        const Real next_root_short( get_root( params.last_short_root ) + params.short_half_wave );

        if( next_root_long < next_root_short )
        {
            next_root_est = next_root_long;

            right_offset = std::min( next_root_short - next_root_est, 
                              params.long_half_wave );

            params.last_long_root = rootList_size() + 1;
        }
        else
        {
            next_root_est = next_root_short;
            
            right_offset = std::min( next_root_long - next_root_est, 
                              params.short_half_wave );

            params.last_short_root = rootList_size() + 1;                        
        }
    }
    
    left_offset = next_root_est - root_n - EPSILON;      

    lower = next_root_est - left_offset;
    upper = next_root_est + safety * right_offset;

    struct root_f_params p = { params.Lm_L, params.h };

    Real f_lower( root_f( lower, &p ) );
    Real f_upper( root_f( upper, &p ) );

    /* set the parity operator for the next root: 
       +1 for an even root.
       -1 for an odd root. */
    const bool parity_op( 2 * ( rootList_size()%2 ) - 1 );

    /* f_lower must have correct sign. */
    assert( f_lower * parity_op < 0 );

    /* If the parity is incorrect -> correct it */
    if( f_upper * parity_op < 0)
    {
        int cntr = 0;

        const Real delta( .1 * std::min(left_offset, right_offset) );
        
        cntr = 0;
        //TODO: find out is upper is left or right from next_root.
        while(f_upper * parity < 0 && cntr++ < 10 )
        {
            upper += delta;
            f_upper = root_f( upper, &p );
        }

        if(cntr >= 10)
        {
            std::cerr << "Failed to straddle root # " << rootList_size() + 1 
                      << " in GF::AbsSinkAbs " << std::endl;
            std::cerr << "f_low(" << lower << ") = " << f_lower << ", f_high(" 
                  << upper << ") = " << f_upper << std::endl;
        }

    }

    return real_pair( lower, upper );
}

/* returns a guess for the number of terms needed for 
   the greensfunction to converge at time t */
uint GreensFunction1DAbsSinkAbs::guess_maxi(Real const& t)
{
    const uint safety(2);

    if (t >= INFINITY)
    {
        return safety;
    }

    const Real D( getD() );

    const Real root0( get_root( 0 ) );
    const Real Dt(D * t);

    const Real thr(exp(- Dt * root0 * root0) * EPSILON * 1e-1);

    if (thr <= 0.0)
    {
        return MAX_TERMS;
    }

    const Real max_root( sqrt(root0 * root0 - log(thr) / Dt) );

    const uint maxi(safety + 
          static_cast<unsigned int>(max_root * ( getLr() + getLl() )  / M_PI));

    return std::min(maxi, MAX_TERMS);
}


/* Standart form of the greensfunction without numerator */
inline Real GreensFunction1DAbsSinkAbs::p_exp_den_i(Real const& t, 
                                                  Real const& root_i, 
                                                  Real const& root_i2) const
{
    return exp( - getD() * root_i2 * t ) / p_denominator_i( root_i );
}


/* Denominator of the greensfunction. */
inline Real GreensFunction1DAbsSinkAbs::p_denominator_i(Real const& root_n) const
{
    const Real Lm( getLr() - getLl() );
    const Real L( getLr() + getLl() );
    
    const Real term1( root_n * L * cos( root_n * L ) + sin( root_n * L ) );
    const Real term2( L * sin( root_n * L ) - Lm * sin( root_n * Lm ) );   

    return getD() * term1 + getk() / 2. * term2;
}


/* Calculates the probability of finding the particle inside the domain
   at time t, the survival probability. */
Real GreensFunction1DAbsSinkAbs::p_survival_old(Real t)
{ 
    const Real D ( getD() );
    const Real k ( getk() ); 
    const Real Lr( getLr() );
    const Real Ll( getLl() );
    const Real L0( getL0() );
    const Real L ( Lr + Ll );


    THROW_UNLESS( std::invalid_argument, t >= 0.0 );

    if (t == 0.0 || (D == 0.0 && v == 0.0) )
    {
	    //particle can't escape.
	    return 1.0;
    }

    Real sum = 0, numerator = 0, prev_term = 0;
    Real term1, term2, term_n, root_n;
    uint n = 0;
    do
    {
        if ( n > MAX_TERMS )
	    {
	        std::cerr << " Too many terms needed for GF1DSink::p_survival. N: "
	              << n << std::endl;
	        break;
	    }
	    prev_term = numerator;
    
        root_n = get_root( n );
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
	        n < MIN_TERMS );

    return 2.0 * sum;
}

/* Calculates survival probability. 
   Switchbox for which greensfunction to use. */
Real GreensFunction3DRadAbs::p_survival(Real t, RealVector& psurvTable) const
{
    THROW_UNLESS( std::invalid_argument, t >= 0.0 );

    Real p;

    const Real D( getD() );
    const Real sigma(getSigma());
    const Real a( geta() );

    const Real L0( getL0() );
    const Real distToa(a - r0);
    const Real distTos(r0 - sigma);

    const Real H(6.0); // a fairly strict criterion for safety.
    const Real maxDist(H * sqrt(2.0 * D * t));

    if (t == 0.0 || (D == 0.0 && v == 0.0) )
    {
	    //particle can't escape.
	    return 1.0;
    }

    /*
    if (distToa > maxDist)
    {
        if (distTos > maxDist) // far from anything; it'll survive.
        {
            if( L0 > maxDist )
                p = 1.0;
        }
        else // close only to s, ignore a
        {
            const Real sigma(this->getSigma());
            const Real kf(this->getkf());
            p = p_survival_irr(t, r0, kf, D, sigma);
        }
    }
    else
    {
        if (distTos > maxDist)  // close only to a.
        {
            p = p_survival_nocollision(t, r0, D, a);
        }
        else  // close to both boundaries.  do the normal calculation.
        {
            const unsigned int maxi(guess_maxi(t));
            
            if (psurvTable.size() < maxi + 1)
            {
                IGNORE_RETURN getAlpha0(maxi);  // this updates the table
                this->createPsurvTable(psurvTable);
            }

            p = funcSum_all(boost::bind(&GreensFunction3DRadAbs::
                                          p_survival_i_exp_table, 
                                          this,
                                          _1, t, psurvTable),
                             maxi);
        }
    }
    */

    const uint maxi( guess_maxi(t) );
            
    if (psurvTable.size() < maxi + 1)
    {
        calculate_n_roots( maxi );  // this updates the table
        createPsurvTable( psurvTable );
    }
    
    p = funcSum_all(boost::bind(&GreensFunction1DAbsSinkAbs::
                                p_survival_i, 
                                this,
                                _1, t, psurvTable),
                    maxi);
    
    return p;
}


/* Calculates the i'th term of the p_survival sum. */
Real p_survival_i( uint const& i, Real const& t, RealVector const& table) const
{
    const Real root_i( get_root( i ) );
    exp( - getD() * t * gsl_pow_2( root_i ) ) * table[ i ];
}


/* Calculates the part of the i'th term of p_surv not dependent on t */
Real p_survival_table_i( Real const& root_i )
{ 
    const Real D ( getD() );
    const Real k ( getk() ); 
    const Real Lr( getLr() );
    const Real Ll( getLl() );
    const Real L0( getL0() );
    const Real L ( Lr + Ll );

    const Real term1( sin( root_n * L ) - sin( root_n * (Lr - L0) ) - sin( root_n * (Ll + L0) ) );
    const Real term2( sin( root_n * Lr ) - sin( root_n * L0 ) - sin( root_n * (Lr - L0) ) );
        
    Real numerator( D * term1 + k * sin( root_n * Ll ) * term2 / root_n ); 
    numerator *= 2;

    return numerator / p_denominator_i( root_i );
}


/* Fills table with terms in the p_survival sum which don't depend on t. */
void GreensFunction1DAbsSinkAbs::createPsurvTable(RealVector& table) const
{
    uint const root_nbr( rootList.size() );
    uint i( table.size() );

    while( i <= root_nbr )
    {
        table.push_back( p_survival_table_i( get_root( i++ ) ) );
    }
}


/* Returns i'th term of prob_r (domain containing r0) function */
Real GreensFunction1DAbsSinkAbs::prob_r_r0_i(uint const& i, Real const& rr, Real const&t)
{
    const Real root_i( get_root(i) );
    const Real LlpL0( Ll + L0 );
    const Real Lrmrr( Lr - rr );
    
    Real numerator( getD() * root_i * sin( root_i * LlpL0 ) + 
                    getk() * sin( root_i * Ll ) * sin( root_i * L0) );

    numerator *= sin( root_i * Lrmrr );
    
    return - 2.0 * p_exp_den_i(t, root_n, gsl_pow_2( root_n ) ) * numerator;
}


/* Returns i'th term of prob_r (domain not containing r0) function */
Real GreensFunction1DAbsSinkAbs::prob_r_nor0_i(uint const& i, Real const& rr, Real const&t)
{
    const Real root_i( get_root(i) );
    const Real LrmL0( Lr - L0 );
    const Real Llprr( Ll + rr );

    const Real numerator( getD() * root_n * sin( root_n * Llprr ) * sin( root_n * LrmL0 ) );

    return -2.0 * p_exp_den_i(t, root_n, gsl_pow_2( root_n ) ) * numerator;
}


/* Calculates the probability density of finding the particle at location r
   at time t. */
Real GreensFunction1DAbsSinkAbs::prob_r(Real r, Real t)
{
    THROW_UNLESS( std::invalid_argument, t >= 0.0 );
    THROW_UNLESS( std::invalid_argument, (r-sigma) >= 0.0 && r <= a && (r0 - sigma) >= 0.0 && r0<=a );
    
    const Real D( getD() );
    const Real Lr( getLr() );
    const Real Ll( getLl() );
    const Real L( Lr + Ll );    
    
    Real L0( getL0() );
    Real rr( r - getrsink() ) ;
    Real p;  

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
    
        p = funcSum_all(boost::bind(&prob_r_r0_i, this, _1, rr, t),
                        maxi);
	}
    else
    {
        if( rr > 0 )
            rr *= -1;

        p = funcSum_all(boost::bind(&prob_r_nor0_i, this, _1, rr, t),
                        maxi);
    }
    
    return p;
}


/* Calculates the probability density of finding the particle at location z at
   timepoint t, given that the particle is still in the domain. */
Real GreensFunction1DAbsSinkAbs::calcpcum(Real r, Real t)
{
    return prob_r(r, t)/p_survival(t);
}


/* Function returns flux at absorbing bounday sigma. */
Real GreensFunction1DAbsSinkAbs::flux_leaves(Real t)
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
Real GreensFunction1DAbsSinkAbs::flux_leavea(Real t)
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
Real GreensFunction1DAbsSinkAbs::flux_tot(Real t)
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
    uint n = 0;
    do
    {
        if ( n > MAX_TERMS )
        {
            std::cerr << " Too many terms needed for GF1DSink::flux_tot. N: "
                      << n << std::endl;
            break;
        }
        prev_term = term_n;
    
        root_n = get_root( n );
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
	        n < MIN_TERMS );

    return  2.0 * D * sum;
}


/* Flux leaving throught absorbing boundary from sub-domain containing r0. */
Real GreensFunction1DAbsSinkAbs::flux_abs_Lr(Real t)
{
    const Real D( getD() );
    const Real k( getk() );    
    const Real Ll( getLl() );
    const Real L0( getL0() );
    const Real LlpL0( Ll + L0 );
    Real root_n;
    
    Real sum = 0, numerator, term_n = 0, prev_term = 0;
    uint n = 0;
    do
    {
        if ( n > MAX_TERMS )
        {
            std::cerr << " Too many terms needed for GF1DSink::flux_abs_Lr. N: "
                      << n << std::endl;
            break;
        }
        prev_term = term_n;
    
        root_n = get_root( n );
        
        numerator = k * sin( root_n * Ll ) * sin ( root_n * L0 ) + D * root_n * sin( root_n * LlpL0 );
        numerator *= root_n;

        term_n = Greens_fn(t, root_n, gsl_pow_2( root_n ) ) * numerator;         
        sum += term_n;

        n++;
    }
    while ( fabs(term_n/sum) > EPSILON  ||
            fabs(prev_term/sum) > EPSILON ||
            n < MIN_TERMS );

    return - D * 2 * sum;
}

/* Flux leaving throught absorbing boundary from sub-domain not containing r0. */
Real GreensFunction1DAbsSinkAbs::flux_abs_Ll(Real t)
{
    const Real D2( gsl_pow_2( getD() ) );
    const Real LrmL0( getLr() - getL0() );

    Real root_n, root_n2;
    
    Real sum = 0, numerator = 0, term_n = 0, prev_term = 0;
    uint n = 0;
    do
    {
        if ( n > MAX_TERMS )
        {
            std::cerr << " Too many terms needed for GF1DSink::flux_abs_Ll. N: "
                      << n << std::endl;
            break;
        }
        prev_term = term_n;
    
        root_n = get_root( n );
        root_n2 = gsl_pow_2( root_n );        
        
        numerator = root_n2 * sin( root_n * LrmL0 );

        term_n = Greens_fn(t, root_n, root_n2 ) * numerator;
        sum += term_n;

        n++;
    }
    while ( fabs(term_n/sum) > EPSILON  ||
            fabs(prev_term/sum) > EPSILON ||
            n < MIN_TERMS );
    
    return 2 * D2 * sum;
}


/* Calculates the probability flux leaving the domain through the sink
   at time t */
Real GreensFunction1DAbsSinkAbs::flux_sink(Real t)
{
    return getk() * prob_r(getrsink(), t);
}


/* Determine which event has occured, an escape or a reaction. Based on the
   fluxes through the boundaries and the sink at the given time. */
GreensFunction1DAbsSinkAbs::EventKind 
GreensFunction1DAbsSinkAbs::drawEventType( Real rnd, Real t )
{
    THROW_UNLESS( std::invalid_argument, rnd < 1.0 && rnd >= 0.0 );
    THROW_UNLESS( std::invalid_argument, t > 0.0 );

    const Real a( geta() );
    const Real sigma( getsigma() );
    const Real r0( getr0() );
    const Real L( a - sigma );

    /* if the sink is impermeable (k==0) or
       the particle is at one the absorbing boundary (sigma or a) => IV_ESCAPE event */
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
Real GreensFunction1DAbsSinkAbs::drawT_f(Real t, void *p)
{
    struct drawT_params *params = (struct drawT_params *)p;
    return params->rnd - params->gf->p_survival( t, params->table );
}


/* Draws the first passage time from the survival probability,
   using an assistance function drawT_f that casts the math. function
   into the form needed by the GSL root solver. */
Real GreensFunction1DAbsSinkAbs::drawTime(Real rnd)
{
    THROW_UNLESS( std::invalid_argument, 0.0 <= rnd && rnd < 1.0 );
  
    const Real a( geta() );    
    const Real r0( getr0() );
    const Real D( getD() );
    const Real Lr( getLr() );
    const Real Ll( getLl() );
    const Real L0( getL0() );
    const Real L( getLr() + getLl() );

    if ( D == 0.0 || L == INFINITY )
    {
	    return INFINITY;
    }

    if ( rnd <= EPSILON || L < 0.0 || fabs( a - r0 ) < EPSILON * L )
    {
	    return 0.0;
    }

    /* the structure to store the numbers to calculate the numbers for 1-S */
    RealVector psurvTable;
    drawT_params params = { this, psurvTable, rnd };

    /* Find a good interval to determine the first passage time.
       First we get the distance to one of the absorbing boundaries or the sink. */
    Real t_guess;
    Real dist( std::min(Lr - L0, Ll + L0) );
    dist = std::min( dist, L0 );

    t_guess = dist * dist / ( 2.0 * D );
    t_guess *= .1;

    const uint n_max_guess( guess_maxi( t_guess ) );
    calculate_n_roots( n_max_guess );

    /* Define the function for the rootfinder */
    gsl_function F;
    F.function = &GreensFunction1DAbsSinkAbs::drawT_f;
    F.params = &parameters;

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

    /* find the intersection on the y-axis between the random number and
       the function */
    
    // define a new solver type brent
    const gsl_root_fsolver_type* solverType( gsl_root_fsolver_brent );
    // make a new solver instance
    gsl_root_fsolver* solver( gsl_root_fsolver_alloc( solverType ) );
    const Real t( findRoot( F, solver, low, high, t_scale*EPSILON, EPSILON,
                            "GreensFunction1DAnsSinkAbs::drawTime" ) );
    // return the drawn time
    return t;
}


/* Retrurns the c.d.f. with respect to the position at time t. */
Real p_int_r(Real const& rr, Real const& t, RealVector& table)
{
    Real p;

    const uint maxi( guess_maxi(t) );       
    if (table.size() < maxi + 1)
    {
        calculate_n_roots( maxi );  // this updates the table
        createP_int_r_Table( t, table );
    }
    
    /* Determine in which part of the domain rr lies, and
       thus which function to use. */
    Real (*p_int_r_i)(uint const&, Real const&, RealVector const& table) = NULL;

    if( rr <= 0 )
        p_int_r_i = &GreensFunction1DAbsSinkAbs::p_int_r_leftdomain;
    else if( rr < obj.getL0() )
            p_int_r_i = &GreensFunction1DAbsSinkAbs::p_int_r_rightdomainA;
        else
            p_int_r_i = &GreensFunction1DAbsSinkAbs::p_int_r_rightdomainB;

    p = funcSum_all(boost::bind(&p_int_r_i, this, _1, rr, t, table),
                    maxi);

    return p
}


void createP_int_r_Table( Real const& t, RealVector& table )
{
    const uint root_nbr( rootList_size() );
    uint i (table.size() );

    while( i <= root_nbr )
    {
        table.push_back( p_exp_den_i(t, get_root( i++ ) ) );
    }
}


/* Function for GFL rootfinder of drawR. */
Real GreensFunction1DAbsSinkAbs::drawR_f(Real rr, void *p)
{
    struct drawR_params *params = (struct drawR_params *)p;
    return params->gf->p_int_r(rr, params->t, params->table) - params->rnd;
}


//Integrated Greens function for rr part of [-Ll, 0]
Real GreensFunction1DAbsSinkAbs::p_int_r_leftdomain(uint const& i, 
                                                    Real const& rr,
                                                    RealVector const& table)
{
    const Real root_i( get_root( i ) );
    const Real LrmL0( getLr() - getL0() );
    const Real Llprr( getLl() + rr );
   
    const Real temp(  2.0 * getD() * sin( root_n * LrmL0 ) * ( cos( root_n * Llprr ) - 1.0 ) );

    return table[i] * temp;
}

//Integrated Greens function for rr part of (0, L0]
Real GreensFunction1DAbsSinkAbs::p_int_r_rightdomainA(uint const& i, 
                                                      Real const& rr, 
                                                      RealVector const& table)
{
    const Real root_i( get_root( i ) );
    const Real LrmL0( getLr() - getL0() );
    const Real Llprr( getLl() + rr );
    const Real root_n_rr( root_n * rr );
    
    const Real temp( getD() * ( cos( root_n * Llprr ) - 1.0 ) + 
                     getk() / root_n * ( cos( root_n_rr ) - 1.0 ) 
                     * sin( root_n * getLl() ) );
    
    return 2.0 * table[i] * sin( root_n * LrmL0 ) * temp;
}

//Integrated Greens function for rr part of (L0, Lr]
Real GreensFunction1DAbsSinkAbs::p_int_r_rightdomainB(uint const& i, 
                                                      Real const& root_n, 
                                                      RealVector const& table)
{
    const Real root_i( get_root( i ) );
    const Real Lr( getLr() );
    const Real Ll( getLl() );
    const Real L0( getL0() );
    const Real L( Lr + Ll );
    const Real LrmL0( Lr - L0 );
    const Real Lrmrr( Lr - rr );
    const Real LlpL0( Ll + L0 );
                   
    const Real term1( sin( root_n * L ) - sin( root_n * LrmL0 ) - 
            sin( root_n * LlpL0 ) * cos( root_n * Lrmrr ) );
            
    const Real term2( sin( root_n * Lr ) - sin( root_n * LrmL0 ) -
            sin( root_n * L0 ) * cos( root_n * Lrmrr ));
        
    const Real temp( getD() * term1 + getk() * sin( root_n * Ll ) * term2 / root_n );

    return 2.0 * table[i] * temp;
}

Real GreensFunction1DAbsSinkAbs::drawR(Real rnd, Real t)
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
    const Real S( p_survival(t) );

    RealVector drawR_table;
    drawR_params parameters = { this, drawR_table, rnd * S };

    const uint n_max_guess( guess_maxi( t ) );

    // define gsl function for rootfinder
    gsl_function F;
    F.function = &drawR_f;
    F.params = &parameters;

    // define a new solver type brent
    const gsl_root_fsolver_type* solverType( gsl_root_fsolver_brent );

    gsl_root_fsolver* solver( gsl_root_fsolver_alloc( solverType ) );
    Real rr( findRoot( F, solver, -Ll, Lr, EPSILON*L, EPSILON,
                            "GreensFunction1DRadAbs::drawR" ) );

    // Convert the position rr to 'world' coordinates and return it.
    if( getr0() < rsink )
        return rsink - rr;
    else
        return rsink + rr;
}

Real GreensFunction1DAbsSinkAbs::drawR_Xn( uint const& i, Real const& t )
{
    const Real root_n( get_root( i ) );
    const Real root_n2( gsl_pow_2( root_n ) );

    return Greens_fn( t, root_n, root_n2 );
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

