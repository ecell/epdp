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
#include "GreensFunction1DRadAbs.hpp"

// This is the appropriate definition of the function defining
// the roots of our Green's functions in GSL.
// Later needed by the rootfinder.
//
// It expects a reaction rate h=k/D already divided by D.
double
GreensFunction1DRadAbs::tan_f (double x, void *p)
{
    // casts the void to the struct pointer
    struct tan_f_params *params = (struct tan_f_params *)p;
    const Real a = (params->a);
    const Real h = (params->h);
    const Real h_a (h*a);

    if ( fabs( h_a ) < 1 )
    {
	// h = k/D
	return 1/tan(x) + (h_a)/x;
    }
    else
    {
	// h = k/D
	return tan(x) + x/(h_a);
    }
    
}

/* Fills the rootList with all the roots of tan(x*a)=-x/h up to n */
void GreensFunction1DRadAbs::calculate_n_roots(uint const& n) const
{
    const Real L( this->geta()-this->getsigma() );
    const Real h( (this->getk()+this->getv()/2.0) / this->getD() );
    // the drift v also comes into this constant, h=(k+v/2)/D
    Real upper, lower, root_i;

    uint i( rootList_size() );

    //No drift, and k = 0, use reflective solution.
    if (h < EPSILON)
    {
        while(i++ < n)
            ad_to_rootList( M_PI * ( n - 1./2 ) / L );
        return;
    }
    
    gsl_function F;
    struct tan_f_params params = { L, h };
     
    F.function = &GreensFunction1DRadAbs::tan_f;
    F.params = &params;

    // define a new solver type brent
    const gsl_root_fsolver_type* solverType( gsl_root_fsolver_brent );
    gsl_root_fsolver* solver( gsl_root_fsolver_alloc( solverType ) );

    /* Find all the roots up to the nth */
    if ( h*L < 1 )
    {
        lower = i*M_PI + 1E-10;
        upper = ( i + 1 ) * M_PI - 1E-10;
    }
    else
    {
        lower = i * M_PI + M_PI_2 + 1E-10;
        upper = ( i + 1 ) * M_PI + M_PI_2 - 1E-10;
    }

    while(i++ < n)
    {
        root_i = findRoot( F, solver, lower, upper, 
                           1.0*EPSILON, EPSILON, "GreensFunction1DRadAbs::root_tan" );

        ad_to_rootList( root_i / L );

        lower += M_PI;
        upper += M_PI;
    }

    gsl_root_fsolver_free( solver );
}


/* returns a guess for the number of terms needed for 
   the greensfunction to converge at time t */
uint GreensFunction1DRadAbs::guess_maxi(Real const& t) const
{
    const uint safety(2);

    if (t >= INFINITY)
    {
        return safety;
    }

    const Real D( getD() );
    const Real L( fabs( geta() - getsigma() ) );

    const Real root0( get_root( 0 ) );
    const Real Dt(D * t);

    const Real thr(exp(- Dt * root0 * root0) * EPSILON * 1e-1);

    if (thr <= 0.0)
    {
        return MAX_TERMS;
    }

    const Real max_root( sqrt(root0 * root0 - log(thr) / Dt) );

    const uint maxi(std::max( safety + 
                              static_cast<uint>
                              (max_root * L  / M_PI),
                              MIN_TERMS )
                    );

    return std::min(maxi, MAX_TERMS);
}


// This is the non-exponential factor in the Green's function sum, not
// including the factor containing the explicit r-dependency (The latter
// is given by the Bn's, see below).
//
// r0 is here still in the interval from 0 to a (and supposed to be the
// starting point of the particle at t0).
//
// The root a_n also must be the specific one for that interval, thus
// the one rescaled by a (see comments in function a_n(n) ).
//
// The factor calculated here is identical for the cases w. or w/o drift,
// only h changes.
Real
GreensFunction1DRadAbs::An (Real root_n) const
{
    const Real h((this->getk()+this->getv()/2.0)/this->getD());
    const Real sigma(this->getsigma());
    const Real L(this->geta()-this->getsigma());
    const Real r0(this->getr0());
    const Real rootn_r0_s = root_n*(r0-sigma);

    return (root_n*cos(rootn_r0_s) + h*sin(rootn_r0_s)) / (h + (root_n*root_n + h*h)*L);
}

// This factor appears in the survival prob.
Real
GreensFunction1DRadAbs::Bn (Real root_n) const
{
    const Real h((this->getk()+this->getv()/2.0)/this->getD());
    const Real k(this->getk());
    const Real D(this->getD());
    const Real v(this->getv());
    const Real sigma(this->getsigma());
    const Real a(this->geta());
    const Real L(this->geta()-this->getsigma());
    
    const Real rootnL(root_n*L);
    const Real rootn2(root_n*root_n);
    const Real h2(h*h);
    const Real v2D(v/2.0/D);

    if(v==0.0) return (h2 - (rootn2 + h2)*cos(rootnL)) / (h*root_n);
    else	return (exp(v2D*sigma)*h*k/D - exp(v2D*a)*(rootn2+h2)*cos(rootnL) ) / (h/root_n*(rootn2+v2D*v2D));
}

// This is the exponential factor in the Green's function sum, also
// appearing in the survival prob. and prop. function.
//
// Also here the root is the one refering to the interval of length L.
Real
GreensFunction1DRadAbs::Cn (Real root_n, Real t)
const
{
    const Real D(this->getD());

    return std::exp(-D*root_n*root_n*t);
}

// Calculates the probability of finding the particle inside the domain
// at time t, the survival probability.
Real
GreensFunction1DRadAbs::p_survival (Real t) const
{
    THROW_UNLESS( std::invalid_argument, t >= 0.0 );
  
    const Real D(this->getD());
    const Real v(this->getv());
    const Real vexpo(-v*v*t/4.0/D - v*r0/2.0/D);
    const Real L(this->geta()-this->getsigma());
    const Real r0_s(this->getr0() - this->getsigma());
    const Real h((this->getk()+this->getv()/2.0)/this->getD());
    const Real r0( getr0() );

    if (t == 0.0 || (D == 0.0 && v == 0.0) )
    {
        // if there was no time or no movement the particle was always
        // in the domain
        return 1.0;
    }

    /* First check if we need full solution. 
       Else we use approximation. */
    const Real distToa( geta() - r0 );
    const Real distTos( r0 - getsigma() );
    const Real maxDist( CUTOFF_H * sqrt(2.0 * D * t) );

    //TODO: No drift included for approximations.
    if( v == 0.0 )
    {
        if( distToa > maxDist ) //Absorbing boundary 'not in sight'.
        {
            if( distTos > maxDist ) //Radiation boundary 'not in sight'.
                return 1.0; //No prob. outflux.
            else
                return XS30(t, distTos, getk(), getD()); //Only radiation BCn.
        }
        else
        {
            if( distTos > maxDist )
                return XS10(t, distToa, getD() ); //Only absorbing BCn.
        }
    }

    Real root_n;
    Real sum = 0, term = 0, term_prev = 0;
    uint n = 0;

    const uint maxi( guess_maxi( t ) );
    calculate_n_roots( maxi );

    if(h != 0.0)
    {
        do
        {
            if ( n >= MAX_TERMS )
            {
                std::cerr << "Too many terms needed for GF1DRad::prob_r. N: "
                          << n << std::endl;
                break;
            }

	        root_n = this->get_root(n);
	        term_prev = term;
    	    term = this->Cn(root_n, t) * this->An(root_n) * this->Bn(root_n);
	        sum += term;
	        n++;
        }
        while ( fabs(term/sum) > EPSILON  ||
                fabs(term_prev/sum) > EPSILON ||
                n < MIN_TERMS);

    }
    else
    {

        do
        {
            if ( n >= MAX_TERMS )
            {
                std::cerr << "Too many terms needed for GF1DRad::prob_r. N: "
                          << n << std::endl;
                break;
            }

	        root_n = this->get_root(n);
	        term_prev = term;
    	    
            term = this->Cn(root_n, t) * cos( root_n * r0_s ) * sin( root_n * L ) / 
                ( L * root_n );

	        sum += term;
	        n++;
        }
        while ( fabs(term/sum) > EPSILON  ||
                fabs(term_prev/sum) > EPSILON ||
                n < MIN_TERMS);

    }
    return 2.0 * exp(vexpo) * sum;
}


// Calculates the probability density of finding the particle at location r
// at time t.
Real
GreensFunction1DRadAbs::prob_r (Real r, Real t)
const
{
    THROW_UNLESS( std::invalid_argument, t >= 0.0 );
    THROW_UNLESS( std::invalid_argument, (r-sigma) >= 0.0 && r <= a 
                  && (r0 - sigma) >= 0.0 && r0<=a );
    
    const Real sigma(this->getsigma());
    const Real a(this->geta());
    const Real L(this->geta()-this->getsigma());
    const Real r0(this->getr0());
    const Real D(this->getD());
    const Real v(this->getv());
    const Real h((this->getk()+this->getv()/2.0)/this->getD());
    
    const Real vexpo(-v*v*t/D/4.0 + v*(r-r0)/D/2.0);

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

    // if r is at the absorbing boundary
    if ( fabs(a-r) < EPSILON*L )
    {
        return 0.0;
    }

    Real root_n, root_n_r_s;
    Real sum = 0, term = 0, prev_term = 0;

    const uint maxi( guess_maxi( t ) );
    calculate_n_roots( maxi );

    uint n = 0;
    do
    {
        if ( n >= MAX_TERMS )
        {
            std::cerr << "Too many terms needed for GF1DRad::prob_r. N: "
                      << n << std::endl;
            break;
        }

        root_n = this->get_root(n);
        root_n_r_s = root_n*(r-sigma);

        prev_term = term;
        term = Cn(root_n, t) * An(root_n) * (h*sin(root_n_r_s) + root_n*cos(root_n_r_s));
        sum += term;

        n++;
    }
    while (fabs(term/sum) > EPSILON*PDENS_TYPICAL || 
           fabs(prev_term/sum) > EPSILON*PDENS_TYPICAL ||
           n < MIN_TERMS );

    return 2.0*exp(vexpo)*sum;
}


// Calculates the probability density of finding the particle at location z at
// timepoint t, given that the particle is still in the domain.
Real
GreensFunction1DRadAbs::calcpcum (Real r, Real t) const
{
    return prob_r(r, t)/p_survival(t);
}


// Calculates the total probability flux leaving the domain at time t
// This is simply the negative of the time derivative of the survival prob.
// at time t [-dS(t')/dt' for t'=t].
Real
GreensFunction1DRadAbs::flux_tot (Real t) const
{
    Real root_n;
    const Real D(this->getD());
    const Real v(this->getv());
    const Real vexpo(-v*v*t/4.0/D - v*r0/2.0/D);

    const Real D2 = D*D;
    const Real v2Dv2D = v*v/4.0/D2;
    double sum = 0, term = 0, prev_term = 0;

    const uint maxi( guess_maxi( t ) );
    calculate_n_roots( maxi );

    uint n = 0;
    do
    {
        if ( n >= MAX_TERMS )
        {
            std::cerr << "Too many terms needed for GF1DRad::flux_tot. N: "
                      << n << std::endl;
            break;
        }

        root_n = this->get_root(n);
        prev_term = term;
        term = (root_n * root_n + v2Dv2D) * Cn(root_n, t) * An(root_n) * Bn(root_n);
        n++;
        sum += term;
    }
    while (fabs(term/sum) > EPSILON*PDENS_TYPICAL ||
           fabs(prev_term/sum) > EPSILON*PDENS_TYPICAL ||
           n < MIN_TERMS );

    return 2.0*D*exp(vexpo)*sum;
}


/* Calculates the probability flux leaving the domain through the radiative
   boundary at time t */
Real
GreensFunction1DRadAbs::flux_rad (Real t) const
{
    return this->getk() * prob_r(this->getsigma(), t);
}


/* Calculates the flux leaving the domain through the radiative boundary as a
   fraction of the total flux. This is the probability that the particle left
   the domain through the radiative boundary instead of the absorbing
   boundary. */
Real
GreensFunction1DRadAbs::fluxRatioRadTot (Real t) const
{
    return flux_rad(t)/flux_tot(t);
}

// Determine which event has occured, an escape or a reaction. Based on the
// fluxes through the boundaries at the given time. Beware: if t is not a
// first passage time you still get an answer!
GreensFunction1DRadAbs::EventKind
GreensFunction1DRadAbs::drawEventType( Real rnd, Real t )
const
{
    THROW_UNLESS( std::invalid_argument, rnd < 1.0 && rnd >= 0.0 );
    THROW_UNLESS( std::invalid_argument, t > 0.0 );
    // if t=0 nothing has happened => no event

    const Real a(this->geta());
    const Real L(this->geta()-this->getsigma());
    const Real r0(this->getr0());

    // if the radiative boundary is impermeable (k==0) or
    // the particle is at the absorbing boundary (at a) => IV_ESCAPE event
    if ( k == 0 || fabs(a-r0) < EPSILON*L )
    {
        return IV_ESCAPE;
    }

    // Else the event is sampled from the flux ratio
    const Real fluxratio (this->fluxRatioRadTot(t));

    if (rnd > fluxratio )
    {
        return IV_ESCAPE;
    }
    else
    {
        return IV_REACTION;
    }
}



/* Returns i'th term of table or fills a table with 
   i exp arguments for drawT. */
Real GreensFunction1DRadAbs::drawT_exponent_table( uint const& i, 
                                                   RealVector& table) const
{
    uint n( table.size() );

    if(i < n)
        return table[i];
    else
    {
        const Real D( getD() );
        const Real v2( gsl_pow_2( getv() ) );
        const Real vexpo_t(-v2/4.0/D);
        Real root_n2;

        while(n++ < i)
        {
            root_n2 = gsl_pow_2( get_root( n ) );
            table.push_back( -D * root_n2 + vexpo_t );
        }
        
        return table.back();
    }
}


/* Returns i'th term of table or fills the table with 
   i numerators for drawT. */
Real GreensFunction1DRadAbs::drawT_Xn_table( uint const& i, 
                                             RealVector& table) const
{
    uint n( table.size() );

    if(i < n)
        return table[i];
    else
    {
        const Real sigma(this->getsigma());
        const Real L(this->geta() - this->getsigma());
        const Real r0(this->getr0());
        const Real D(this->getD());
        const Real v(this->getv());
        const Real h((getk()+getv()/2.0)/getD());

        Real root_n, root_n2, root_n_r0_s, root_n_L, h_root_n;	  

        if( v == 0.0)
        {
            while(n++ < i)
            {
                root_n      = get_root(n);	
                root_n2	    = root_n * root_n;
                root_n_r0_s = root_n * (r0-sigma);
                root_n_L    = root_n * L;
                h_root_n    = h / root_n;

                table.push_back( (h*sin(root_n_r0_s) + root_n*cos(root_n_r0_s)) 
                                 / (L*(root_n2+h*h)+h) * ( h_root_n + sin(root_n_L) 
                                                           - h_root_n*cos(root_n_L) ) );
            }
        }
        else
        {
            const Real v2D(v/2.0/D);
            const Real exp_av2D(exp(a*v2D));
            const Real exp_sigmav2D(exp(sigma*v2D));

            while(n++ < i)
            {
                root_n      = get_root(n);	
                root_n2	    = root_n * root_n;
                root_n_r0_s = root_n * (r0-sigma);
                root_n_L    = root_n * L;
                h_root_n    = h / root_n;

                table.push_back( (h*sin(root_n_r0_s) + root_n*cos(root_n_r0_s)) / 
                                 (L*(root_n2+h*h)+h) * (exp_sigmav2D*h*k/D - 
                                                        exp_av2D*(root_n2+h*h)*
                                                        cos(root_n_L)) / 
                                 (h_root_n * (root_n2 + v2D*v2D))  
                                 );
            }
        }
        return table.back();
    }
}


// This function is needed to cast the math. form of the function
// into the form needed by the GSL root solver.
double GreensFunction1DRadAbs::drawT_f (double t, void *p)
{
    // casts p to type 'struct drawT_params *'
    struct drawT_params *params = (struct drawT_params *)p;
    GreensFunction1DRadAbs const* gf = params->gf;
    const Real r0( gf->getr0() );

    /* First check if we need full solution. 
       Else we use approximation. */
    const Real distToa( gf->geta() - r0 );
    const Real distTos( r0 - gf->getsigma() );
    const Real maxDist( CUTOFF_H * sqrt(2.0 * gf->getD() * t) );
    
    //TODO: No drift included for approximations.
    if( gf->getv() == 0.0 )
    {
        if( distToa > maxDist ) //Absorbing boundary 'not in sight'.
        {
            if( distTos > maxDist )//And radiation boundary 'not in sight'.
                return params->rnd - 1.0; //No prob. outflux.
            else
                return params->rnd - XS30(t, distTos, gf->getk(), gf->getD()); //Only radiation BCn.
        }
        else
        {
            if( distTos > maxDist )
                return params->rnd - XS10(t, distToa, gf->getD() ); //Only absorbing BCn.
        }
    }


    Real Xn = 0, exponent = 0;
    Real prefactor = params->prefactor;

    const uint maxi( gf->guess_maxi( t ) );
    
    if( params->exponent_table.size() < maxi + 1 &&
        params->Xn_table.size() < maxi + 1 )
    {
        params->gf->calculate_n_roots( maxi );
        IGNORE_RETURN gf->drawT_exponent_table( maxi, params->exponent_table );
        IGNORE_RETURN gf->drawT_Xn_table( maxi, params->Xn_table );
    }

    Real sum = 0, term = 0, prev_term = 0;
    uint n = 0;
    do
    {
        if ( n >= maxi )
        {
            std::cerr << "Too many terms needed for GF1DRad::DrawTime. N: "
                      << n << ", conv arg: " << exponent * t <<  std::endl;
            break;
        }
        prev_term = term;
        
        Xn = params->gf->drawT_Xn_table(n, params->Xn_table);
        exponent = params->gf->drawT_exponent_table(n, params->exponent_table);
        term = Xn * exp(exponent * t);

        sum += term;
        n++;
    }
    while (fabs(term/sum) > EPSILON*1.0 ||
           fabs(prev_term/sum) > EPSILON*1.0 ||
           n < MIN_TERMS );

    // find the intersection with the random number
    return params->rnd - prefactor*sum;
}

// Draws the first passage time from the survival probability,
// using an assistance function drawT_f that casts the math. function
// into the form needed by the GSL root solver.
Real
GreensFunction1DRadAbs::drawTime (Real rnd) const
{
    THROW_UNLESS( std::invalid_argument, 0.0 <= rnd && rnd < 1.0 );
  
    const Real sigma(this->getsigma());
    const Real a(this->geta());
    const Real L(this->geta()-this->getsigma());
    const Real r0(this->getr0());
    const Real k(this->getk());
    const Real D(this->getD());
    const Real v(this->getv());   

    if ( D == 0.0 || L == INFINITY )
    {
        return INFINITY;
    }

    if ( rnd > 1 - EPSILON || L < 0.0 || fabs(a-r0) < EPSILON*L )
    {
        return 0.0;
    }

    // some temporary variables
    const Real vexpo_pref(-v*r0/2.0/D);
    Real prefactor;

    /* Find a good interval to determine the first passage time. */
    Real t_guess;
    Real dist;

    if (k != 0)
    {
        dist = std::min(a - r0, r0 - sigma);
    }
    else
    {
        dist = a - r0;
    }

    t_guess = dist * dist / ( 2.0 * D );
    t_guess *= .1;

    // A different guess has to be made in case of nonzero drift to account for the displacement due to it
    // TODO: This does not work properly in this case yet, but we don't know why...
    // When drifting towards the closest boundary
    //if( (r0 >= a/2.0 && v > 0.0) || (r0 <= a/2.0 && v < 0.0) )	t_guess = sqrt(D*D/(v*v*v*v)+dist*dist/(v*v)) - D/(v*v);
    // When drifting away from the closest boundary
    //if( ( r0 < a/2.0 && v > 0.0) || ( r0 > a/2.0 && v < 0.0) )	t_guess = D/(v*v) - sqrt(D*D/(v*v*v*v)-dist*dist/(v*v));  

    /* Set params structure. */
    RealVector exponent_table;
    RealVector Xn_table;
    struct drawT_params parameters = {this, exponent_table, Xn_table};
    
    // the prefactor of the sum is also different in case of drift<>0 :
    if(v == 0)
        prefactor = 2.0;
    else
        prefactor = 2.0 * exp(vexpo_pref);

    parameters.prefactor  = prefactor;
    
    // store the random number for the probability
    parameters.rnd = rnd;
    // store the number of terms used
    parameters.tscale = this->t_scale;

    // Define the function for the rootfinder
    gsl_function F;
    F.function = &GreensFunction1DRadAbs::drawT_f;
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
                std::cerr << "GF1DRad::drawTime Couldn't adjust high. F("
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
                std::cerr << "GF1DRad::drawTime Couldn't adjust low. F(" << low << ") = "
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
                            "GreensFunction1DRadAbs::drawTime" ) );

    // return the drawn time
    return t;
}

double
GreensFunction1DRadAbs::drawR_f (double z, void *p)
{
    // casts p to type 'struct drawR_params *'
    struct drawR_params *params = (struct drawR_params *)p;
    Real v2D 		= params->H[0];	// = v2D = v/(2D)
    Real costerm 	= params->H[1];	// = k/D
    Real sinterm 	= params->H[2];	// = h*v2D
    Real sigma 		= params->H[3];	// = sigma
    uint  terms = params->terms;

    Real expsigma(exp(sigma*v2D));
    Real zs(z-sigma);
    
    Real sum = 0, term = 0, prev_term = 0;
    Real root_n, S_Cn_root_n;
    
    uint n = 0;
    do
    {
        if ( n >= terms )
        {
            std::cerr << "GF1DRad: Too many terms needed for DrawR. N: "
                      << n << std::endl;
            break;
        }
        prev_term = term;
        
        S_Cn_root_n = params->S_Cn_root_n[n];
        root_n  = params->root_n[n];
        term = S_Cn_root_n * ( expsigma*costerm - exp(v2D*z)*
                               ( costerm*cos(root_n*zs) - 
                                 (root_n+sinterm/root_n)*sin(root_n*zs) ));

        sum += term;
        n++;
    }
    while (fabs(term/sum) > EPSILON*1.0 ||
           fabs(prev_term/sum) > EPSILON*1.0 ||
           n < MIN_TERMS );

    // Find the intersection with the random number
    return sum - params->rnd;
}

Real
GreensFunction1DRadAbs::drawR (Real rnd, Real t) const
{
    THROW_UNLESS( std::invalid_argument, 0.0 <= rnd && rnd < 1.0 );
    THROW_UNLESS( std::invalid_argument, t >= 0.0 );
    
    const Real sigma(this->getsigma());
    const Real a(this->geta());
    const Real L(this->geta()-this->getsigma());
    const Real r0(this->getr0());
    const Real D(this->getD());
    const Real v(this->getv());
    const Real k(this->getk());
    const Real h((this->getk()+this->getv()/2.0)/this->getD());
    

    if (t==0.0 || (D==0.0 && v==0.0) )
    {
        // the trivial case
        //return r0*this->l_scale;	// renormalized version, discontinued
        return r0;
    }
    if ( a<0.0 )
    {
        // if the domain had zero size
        return 0.0;
    }

    // the structure to store the numbers to calculate the numbers for 1-S
    struct drawR_params parameters;
    double root_n = 0;
    double S_Cn_root_n;
    double root_n2, root_n_r0_s;

    const Real vexpo(-v*v*t/4.0/D - v*r0/2.0/D); // exponent of the drift-prefactor, same as in survival prob.
    const Real v2D(v/2.0/D);
    const Real v2Dv2D(v2D*v2D);
    const Real S = 2.0*exp(vexpo)/p_survival(t); 

    assert(p_survival(t) >= 0.0);

    const uint maxi( guess_maxi( t ) );
    calculate_n_roots( maxi );

    // produce the coefficients and the terms in the exponent and put them
    // in the params structure
    for (uint n = 0; n < maxi; n++)
    {
	    root_n = this->get_root(n);  // get the n-th root of tan(alfa*a)=alfa/-k
       	root_n2 = root_n * root_n;
       	root_n_r0_s = root_n * (r0-sigma);
       	S_Cn_root_n =	S * exp(-D*root_n2*t)
	         * (root_n*cos(root_n_r0_s) + h*sin(root_n_r0_s)) / (L*(root_n2 + h*h) + h)
	         * root_n / (root_n2 + v2Dv2D);

      	// store the coefficients in the structure
	    parameters.root_n[n] = root_n;
        // also store the values for the exponent
        parameters.S_Cn_root_n[n] = S_Cn_root_n;
    }
    
    // also store constant prefactors that appear in the calculation of the
    // r-dependent terms
    parameters.H[0] = v2D;		// appears together with z in one of the prefactors
    parameters.H[1] = k/D;		// further constant terms of the cosine prefactor
    parameters.H[2] = h*v2D;		// further constant terms of the sine prefactor
    parameters.H[3] = sigma;

    // define gsl function for rootfinder
    gsl_function F;
    F.function = &GreensFunction1DRadAbs::drawR_f;
    F.params = &parameters;

    // store the random number for the probability
    parameters.rnd = rnd;
    // store the number of terms used
    parameters.terms = maxi;

    // define a new solver type brent
    const gsl_root_fsolver_type* solverType( gsl_root_fsolver_brent );
    gsl_root_fsolver* solver( gsl_root_fsolver_alloc( solverType ) );

    Real r( findRoot( F, solver, sigma, a, EPSILON*L, EPSILON,
                            "GreensFunction1DRadAbs::drawR" ) );

    // return the drawn position
    return r;
}

std::string GreensFunction1DRadAbs::dump() const
{
    std::ostringstream ss;
    ss << "D = " << this->getD() << ", sigma = " << this->getsigma() <<
        ", a = " << this->geta() <<
        ", r0 = " << this->getr0() <<
        ", k = " << this->getk() << std::endl;
    return ss.str();
}

