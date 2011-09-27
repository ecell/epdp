#if !defined( __FIRSTPASSAGEPAIRGREENSFUNCTION2D_HPP )
#define __FIRSTPASSAGEPAIRGREENSFUNCTION2D_HPP 

#include <vector>
#include <boost/tuple/tuple.hpp>
#include <boost/function.hpp>
#include <boost/array.hpp>

#include <gsl/gsl_roots.h>

#include "PairGreensFunction.hpp"


class GreensFunction2DRadAbs
    :
    public PairGreensFunction
{

public:
    // defines vector from template defined in standard template library,
    // gives it's membervalues type Real. std:: clarifies the std namespace.
    typedef std::vector<Real> RealVector;

private:
    // Error tolerance used by default.
    static const Real TOLERANCE = 1e-8;

    // SphericalBesselGenerator's accuracy, used by some
    // theta-related calculations.

    static const Real MIN_T_FACTOR = 1e-8;

    static const Real L_TYPICAL = 1E-7;
    static const Real T_TYPICAL = 1E-5;
    static const Real EPSILON = 1E-12;

    static const unsigned int MAX_ORDER = 30;		// The maximum number of m terms
    static const unsigned int MAX_ALPHA_SEQ = 2000;	// The maximum number of n terms

public:

    const char* getName() const
    {
        return "GreensFunction2DRadAbs";
    }
    
    GreensFunction2DRadAbs( const Real D, 
				    const Real kf,
               			    const Real r0, 
				    const Real Sigma,
		                    const Real a );
  
    virtual ~GreensFunction2DRadAbs();

    const Real geth() const
    {
	return this->h;
    }

    const Real geta() const
    {
	return this->a;
    }

    virtual Real drawTime( const Real rnd) const;

    virtual EventKind drawEventType( const Real rnd, 
				   const Real t ) const;
    
    virtual Real drawR( const Real rnd, 
		      const Real t ) const;
    
    virtual Real drawTheta( const Real rnd,
			  const Real r, 
			  const Real t ) const;
    
    
    const Real f_alpha0( const Real alpha ) const;
  
    const Real f_alpha( const Real alpha, const Integer n ) const;

    
    const Real p_survival( const Real t) const;

    const Real p_survival_table( const Real t,
				 RealVector& table ) const;


    const Real leaves( const Real t) const;

    const Real leavea( const Real t) const;

    const Real p_m( const Integer n, const Real r, const Real t ) const;

    const Real dp_m_at_a( const Integer m, const Real t ) const;


    const Real p_m_alpha( const unsigned int n,
			  const unsigned int m,
			  const Real r, 
			  const Real t ) const;

    const Real dp_m_alpha_at_a( const unsigned int n,
				const unsigned int m,
				const Real t ) const;

    // methods below are kept public for debugging purpose.

    std::string dump() const;

    const Real alphaOffset( const unsigned int n ) const;

    const Real alpha0_i( const Real previous ) const;

    const Real alpha_i( const Real offset, const Integer n ) const;

    const Real p_survival_i( const Real alpha) const;

    const Real calc_A_i_0( const Real alpha) const;

    const Real leaves_i( const Real alpha) const;

    const boost::tuple<Real,Real,Real> Y0J0J1_constants ( const Real alpha,
                                                          const Real t) const;

    const Real getAlpha( const size_t n, const RealVector::size_type i ) const;

    const Real getAlpha0( const RealVector::size_type i ) const;

protected:

    void clearAlphaTable() const;


    RealVector& getAlphaTable( const size_t n ) const
    {
	return this->alphaTable[n];
    }

    const Real p_int_r_table( const Real r,
				const RealVector& Y0_aAnTable,
				const RealVector& J0_aAnTable,
				const RealVector& Y0J1J0Y1Table ) const;

    const Real ip_theta_table( const Real theta,
			       const RealVector& p_nTable ) const;

    const Real p_survival_i_exp_table( const unsigned int i,
				       const Real t,
				       const RealVector& table ) const;

    const Real leavea_i_exp( const unsigned int i,
			     const Real alpha) const;

    const Real leaves_i_exp( const unsigned int i,
			     const Real alpha) const;

    const Real ip_theta_n( const unsigned int m,
			   const Real theta,
			   const RealVector& p_nTable ) const;


    const Real p_int_r_i_exp_table( const unsigned int i,
					const Real r,
					const RealVector& Y0_aAnTable,
					const RealVector& J0_aAnTable,
					const RealVector& Y0J1J0Y1Table ) const;

    void createPsurvTable( RealVector& table) const; 

    void createY0J0Tables( RealVector& Y0_Table, RealVector& J0_Table, RealVector& Y0J1J0Y1_Table,
				const Real t ) const;

    void makep_mTable( RealVector& p_mTable,
		       const Real r, 
		       const Real t ) const;
    
    void makedp_m_at_aTable( RealVector& p_mTable,
			     const Real t ) const;

    const unsigned int guess_maxi( const Real t ) const;

    struct f_alpha0_aux_params
    { 
	const GreensFunction2DRadAbs* const gf;
	const Real value;
    };

    static const Real 
    f_alpha0_aux_F( const Real alpha,
		    const f_alpha0_aux_params* const params );


    struct f_alpha_aux_params
    { 
	const GreensFunction2DRadAbs* const gf;
	const Integer n;
	Real value;
    };

    static const Real 
    f_alpha_aux_F( const Real alpha,
		   const f_alpha_aux_params* const params );

    struct p_survival_table_params
    { 
	const GreensFunction2DRadAbs* const gf;
//	const Real r0;
	RealVector& table;
	const Real rnd;
    };

    static const Real 
    p_survival_table_F( const Real t,
                        const p_survival_table_params* const params );

    struct p_int_r_params
    { 
	const GreensFunction2DRadAbs* const gf;
	const Real t;
//	const Real r0;
	const RealVector& Y0_aAnTable;
	const RealVector& J0_aAnTable;
	const RealVector& Y0J1J0Y1Table;
	const Real rnd;
    };

    static const Real 
    p_int_r_F( const Real r,
	       const p_int_r_params* const params );

    struct ip_theta_params
    { 
	const GreensFunction2DRadAbs* const gf;
	const Real r;
//	const Real r0;
	const Real t;
	const RealVector& p_nTable;
	const Real value;
    };

    static const Real 
    ip_theta_F( const Real theta,
		const ip_theta_params* const params );

private:
    
    const Real h;

    mutable boost::array<Real,MAX_ORDER+1> alphaOffsetTable;
    mutable boost::array<RealVector,MAX_ORDER+1> alphaTable;

    Real a;

};



#endif // __FIRSTPASSAGEPAIRGREENSFUNCTION2D_HPP
