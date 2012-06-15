#ifndef STRUCTFUNC_CPP
#define STRUCTFUNC_CPP

//#include "StructureFunctions.hpp" // TODO Why the hell does this not work, man???
#include "CuboidalRegion.hpp"
#include "SphericalSurface.hpp"
#include "CylindricalSurface.hpp"
#include "DiskSurface.hpp"
#include "PlanarSurface.hpp"
#include "exceptions.hpp"


/****************************/
/*** INTERACTION HANDLERS ***/
/****************************/

/******************************/
/* Coming from CuboidalRegion */
/******************************/
template <typename Ttraits_>
inline std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>
get_pos_sid_pair( CuboidalRegion<Ttraits_>              const& origin_structure,
                  CuboidalRegion<Ttraits_>              const& target_structure,
                  typename Ttraits_::position_type      const& old_pos          )
{
    typedef typename Ttraits_::structure_id_type        structure_id_type;
    typedef typename Ttraits_::position_type            position_type;

    /*** COMBINATION NOT SUPPORTED ***/
    throw illegal_propagation_attempt("Interaction between combination of origin structure and target structure not supported.");
    
    return std::make_pair(position_type(), structure_id_type());
};

template <typename Ttraits_>
inline std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>
get_pos_sid_pair( CuboidalRegion<Ttraits_>              const& origin_structure,
                  SphericalSurface<Ttraits_>            const& target_structure,
                  typename Ttraits_::position_type      const& old_pos          )
{
    typedef typename Ttraits_::structure_id_type        structure_id_type;
    typedef typename Ttraits_::position_type            position_type;

    /*** COMBINATION NOT SUPPORTED ***/
    throw illegal_propagation_attempt("Interaction between combination of origin structure and target structure not supported.");
    
    return std::make_pair(position_type(), structure_id_type());
};

template <typename Ttraits_>
inline std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>
get_pos_sid_pair( CuboidalRegion<Ttraits_>              const& origin_structure,
                  CylindricalSurface<Ttraits_>          const& target_structure,
                  typename Ttraits_::position_type      const& old_pos          )
{
    typedef typename Ttraits_::structure_id_type        structure_id_type;
    typedef typename Ttraits_::position_type            position_type;

    /*** COMBINATION NOT SUPPORTED ***/
    throw illegal_propagation_attempt("Interaction between combination of origin structure and target structure not supported.");
    
    return std::make_pair(position_type(), structure_id_type());
};

template <typename Ttraits_>
inline std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>
get_pos_sid_pair( CuboidalRegion<Ttraits_>              const& origin_structure,
                  DiskSurface<Ttraits_>                 const& target_structure,
                  typename Ttraits_::position_type      const& old_pos          )
{
    typedef typename Ttraits_::structure_id_type        structure_id_type;
    typedef typename Ttraits_::position_type            position_type;

    /*** COMBINATION NOT SUPPORTED ***/
    throw illegal_propagation_attempt("Interaction between combination of origin structure and target structure not supported.");
    
    return std::make_pair(position_type(), structure_id_type());
};

template <typename Ttraits_>
inline std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>
get_pos_sid_pair( CuboidalRegion<Ttraits_>              const& origin_structure,
                  PlanarSurface<Ttraits_>               const& target_structure,
                  typename Ttraits_::position_type      const& old_pos          )
{
    typedef typename Ttraits_::structure_id_type        structure_id_type;
    typedef typename Ttraits_::position_type            position_type;

    /*** COMBINATION NOT SUPPORTED ***/
    throw illegal_propagation_attempt("Interaction between combination of origin structure and target structure not supported.");
    
    return std::make_pair(position_type(), structure_id_type());
};

/********************************/
/* Coming from SphericalSurface */
/********************************/
template <typename Ttraits_>
inline std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>
get_pos_sid_pair( SphericalSurface<Ttraits_>            const& origin_structure,
                  CuboidalRegion<Ttraits_>              const& target_structure,
                  typename Ttraits_::position_type      const& old_pos          )
{
    typedef typename Ttraits_::structure_id_type        structure_id_type;
    typedef typename Ttraits_::position_type            position_type;

    /*** COMBINATION NOT SUPPORTED ***/
    throw illegal_propagation_attempt("Interaction between combination of origin structure and target structure not supported.");
    
    return std::make_pair(position_type(), structure_id_type());
};

template <typename Ttraits_>
inline std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>
get_pos_sid_pair( SphericalSurface<Ttraits_>            const& origin_structure,
                  SphericalSurface<Ttraits_>            const& target_structure,
                  typename Ttraits_::position_type      const& old_pos          )
{
    typedef typename Ttraits_::structure_id_type        structure_id_type;
    typedef typename Ttraits_::position_type            position_type;

    /*** COMBINATION NOT SUPPORTED ***/
    throw illegal_propagation_attempt("Interaction between combination of origin structure and target structure not supported.");
    
    return std::make_pair(position_type(), structure_id_type());
};

template <typename Ttraits_>
inline std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>
get_pos_sid_pair( SphericalSurface<Ttraits_>            const& origin_structure,
                  CylindricalSurface<Ttraits_>          const& target_structure,
                  typename Ttraits_::position_type      const& old_pos          )
{
    typedef typename Ttraits_::structure_id_type        structure_id_type;
    typedef typename Ttraits_::position_type            position_type;

    /*** COMBINATION NOT SUPPORTED ***/
    throw illegal_propagation_attempt("Interaction between combination of origin structure and target structure not supported.");
    
    return std::make_pair(position_type(), structure_id_type());
};

template <typename Ttraits_>
inline std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>
get_pos_sid_pair( SphericalSurface<Ttraits_>            const& origin_structure,
                  DiskSurface<Ttraits_>                 const& target_structure,
                  typename Ttraits_::position_type      const& old_pos          )
{
    typedef typename Ttraits_::structure_id_type        structure_id_type;
    typedef typename Ttraits_::position_type            position_type;

    /*** COMBINATION NOT SUPPORTED ***/
    throw illegal_propagation_attempt("Interaction between combination of origin structure and target structure not supported.");
    
    return std::make_pair(position_type(), structure_id_type());
};

template <typename Ttraits_>
inline std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>
get_pos_sid_pair( SphericalSurface<Ttraits_>            const& origin_structure,
                  PlanarSurface<Ttraits_>               const& target_structure,
                  typename Ttraits_::position_type      const& old_pos          )
{
    typedef typename Ttraits_::structure_id_type        structure_id_type;
    typedef typename Ttraits_::position_type            position_type;

    /*** COMBINATION NOT SUPPORTED ***/
    throw illegal_propagation_attempt("Interaction between combination of origin structure and target structure not supported.");
    
    return std::make_pair(position_type(), structure_id_type());
};

/**********************************/
/* Coming from CylindricalSurface */
/**********************************/
template <typename Ttraits_>
inline std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>
get_pos_sid_pair( CylindricalSurface<Ttraits_>          const& origin_structure,
                  CuboidalRegion<Ttraits_>              const& target_structure,
                  typename Ttraits_::position_type      const& old_pos          )
{
    typedef typename Ttraits_::structure_id_type        structure_id_type;
    typedef typename Ttraits_::position_type            position_type;

    /*** COMBINATION NOT SUPPORTED ***/
    throw illegal_propagation_attempt("Interaction between combination of origin structure and target structure not supported.");
    
    return std::make_pair(position_type(), structure_id_type());
};

template <typename Ttraits_>
inline std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>
get_pos_sid_pair( CylindricalSurface<Ttraits_>          const& origin_structure,
                  CylindricalSurface<Ttraits_>          const& target_structure,
                  typename Ttraits_::position_type      const& old_pos          )
{
    typedef typename Ttraits_::structure_id_type        structure_id_type;
    typedef typename Ttraits_::position_type            position_type;

    /*** COMBINATION NOT SUPPORTED ***/
    throw illegal_propagation_attempt("Interaction between combination of origin structure and target structure not supported.");
    
    return std::make_pair(position_type(), structure_id_type());
};

template <typename Ttraits_>
inline std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>
get_pos_sid_pair( CylindricalSurface<Ttraits_>          const& origin_structure,
                  SphericalSurface<Ttraits_>            const& target_structure,
                  typename Ttraits_::position_type      const& old_pos          )
{
    typedef typename Ttraits_::structure_id_type        structure_id_type;
    typedef typename Ttraits_::position_type            position_type;

    /*** COMBINATION NOT SUPPORTED ***/
    throw illegal_propagation_attempt("Interaction between combination of origin structure and target structure not supported.");
    
    return std::make_pair(position_type(), structure_id_type());
};

template <typename Ttraits_>
inline std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>
get_pos_sid_pair( CylindricalSurface<Ttraits_>          const& origin_structure,
                  DiskSurface<Ttraits_>                 const& target_structure,
                  typename Ttraits_::position_type      const& old_pos          )
{
    typedef typename Ttraits_::structure_id_type        structure_id_type;
    typedef typename Ttraits_::position_type            position_type;

    /*** COMBINATION NOT SUPPORTED ***/
    throw illegal_propagation_attempt("Interaction between combination of origin structure and target structure not supported.");
    
    return std::make_pair(position_type(), structure_id_type());
};

template <typename Ttraits_>
inline std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>
get_pos_sid_pair( CylindricalSurface<Ttraits_>          const& origin_structure,
                  PlanarSurface<Ttraits_>               const& target_structure,
                  typename Ttraits_::position_type      const& old_pos          )
{
    typedef typename Ttraits_::structure_id_type        structure_id_type;
    typedef typename Ttraits_::position_type            position_type;

    /*** COMBINATION NOT SUPPORTED ***/
    throw illegal_propagation_attempt("Interaction between combination of origin structure and target structure not supported.");
    
    return std::make_pair(position_type(), structure_id_type());
};

/***************************/
/* Coming from DiskSurface */
/***************************/
template <typename Ttraits_>
inline std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>
get_pos_sid_pair( DiskSurface<Ttraits_>                 const& origin_structure,
                  CuboidalRegion<Ttraits_>              const& target_structure,
                  typename Ttraits_::position_type      const& old_pos          )
{
    typedef typename Ttraits_::structure_id_type        structure_id_type;
    typedef typename Ttraits_::position_type            position_type;

    /*** COMBINATION NOT SUPPORTED ***/
    throw illegal_propagation_attempt("Interaction between combination of origin structure and target structure not supported.");
    
    return std::make_pair(position_type(), structure_id_type());
};

template <typename Ttraits_>
inline std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>
get_pos_sid_pair( DiskSurface<Ttraits_>                 const& origin_structure,
                  SphericalSurface<Ttraits_>            const& target_structure,
                  typename Ttraits_::position_type      const& old_pos          )
{
    typedef typename Ttraits_::structure_id_type        structure_id_type;
    typedef typename Ttraits_::position_type            position_type;

    /*** COMBINATION NOT SUPPORTED ***/
    throw illegal_propagation_attempt("Interaction between combination of origin structure and target structure not supported.");
    
    return std::make_pair(position_type(), structure_id_type());
};

template <typename Ttraits_>
inline std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>
get_pos_sid_pair( DiskSurface<Ttraits_>                 const& origin_structure,
                  CylindricalSurface<Ttraits_>          const& target_structure,
                  typename Ttraits_::position_type      const& old_pos          )
{
    typedef typename Ttraits_::structure_id_type        structure_id_type;
    typedef typename Ttraits_::position_type            position_type;

    /*** COMBINATION NOT SUPPORTED ***/
    throw illegal_propagation_attempt("Interaction between combination of origin structure and target structure not supported.");
    
    return std::make_pair(position_type(), structure_id_type());
};

template <typename Ttraits_>
inline std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>
get_pos_sid_pair( DiskSurface<Ttraits_>                 const& origin_structure,
                  DiskSurface<Ttraits_>                 const& target_structure,
                  typename Ttraits_::position_type      const& old_pos          )
{
    typedef typename Ttraits_::structure_id_type        structure_id_type;
    typedef typename Ttraits_::position_type            position_type;

    /*** COMBINATION NOT SUPPORTED ***/
    throw illegal_propagation_attempt("Interaction between combination of origin structure and target structure not supported.");
    
    return std::make_pair(position_type(), structure_id_type());
};

template <typename Ttraits_>
inline std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>
get_pos_sid_pair( DiskSurface<Ttraits_>                 const& origin_structure,
                  PlanarSurface<Ttraits_>               const& target_structure,
                  typename Ttraits_::position_type      const& old_pos          )
{
    typedef typename Ttraits_::structure_id_type        structure_id_type;
    typedef typename Ttraits_::position_type            position_type;

    /*** COMBINATION NOT SUPPORTED ***/
    throw illegal_propagation_attempt("Interaction between combination of origin structure and target structure not supported.");
    
    return std::make_pair(position_type(), structure_id_type());
};

/*****************************/
/* Coming from PlanarSurface */
/*****************************/
template <typename Ttraits_>
inline std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>
get_pos_sid_pair( PlanarSurface<Ttraits_>               const& origin_structure,
                  CuboidalRegion<Ttraits_>              const& target_structure,
                  typename Ttraits_::position_type      const& old_pos          )
{
    typedef typename Ttraits_::structure_id_type        structure_id_type;
    typedef typename Ttraits_::position_type            position_type;

    /*** COMBINATION NOT SUPPORTED ***/
    throw illegal_propagation_attempt("Interaction between combination of origin structure and target structure not supported.");
    
    return std::make_pair(position_type(), structure_id_type());
};

template <typename Ttraits_>
inline std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>
get_pos_sid_pair( PlanarSurface<Ttraits_>               const& origin_structure,
                  SphericalSurface<Ttraits_>            const& target_structure,
                  typename Ttraits_::position_type      const& old_pos          )
{
    typedef typename Ttraits_::structure_id_type        structure_id_type;
    typedef typename Ttraits_::position_type            position_type;

    /*** COMBINATION NOT SUPPORTED ***/
    throw illegal_propagation_attempt("Interaction between combination of origin structure and target structure not supported.");
    
    return std::make_pair(position_type(), structure_id_type());
};

template <typename Ttraits_>
inline std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>
get_pos_sid_pair( PlanarSurface<Ttraits_>               const& origin_structure,
                  CylindricalSurface<Ttraits_>          const& target_structure,
                  typename Ttraits_::position_type      const& old_pos          )
{
    typedef typename Ttraits_::structure_id_type        structure_id_type;
    typedef typename Ttraits_::position_type            position_type;

    /*** COMBINATION NOT SUPPORTED ***/
    throw illegal_propagation_attempt("Interaction between combination of origin structure and target structure not supported.");
    
    return std::make_pair(position_type(), structure_id_type());
};

template <typename Ttraits_>
inline std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>
get_pos_sid_pair( PlanarSurface<Ttraits_>               const& origin_structure,
                  DiskSurface<Ttraits_>                 const& target_structure,
                  typename Ttraits_::position_type      const& old_pos          )
{
    typedef typename Ttraits_::structure_id_type        structure_id_type;
    typedef typename Ttraits_::position_type            position_type;

    /*** COMBINATION NOT SUPPORTED ***/
    throw illegal_propagation_attempt("Interaction between combination of origin structure and target structure not supported.");
    
    return std::make_pair(position_type(), structure_id_type());
};

template <typename Ttraits_>
inline std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>
get_pos_sid_pair( PlanarSurface<Ttraits_>               const& origin_structure,
                  PlanarSurface<Ttraits_>               const& target_structure,
                  typename Ttraits_::position_type      const& old_pos          )
{
    typedef typename Ttraits_::structure_id_type        structure_id_type;
    typedef typename Ttraits_::position_type            position_type;

    /*** COMBINATION NOT SUPPORTED ***/
    throw illegal_propagation_attempt("Interaction between combination of origin structure and target structure not supported.");
    
    return std::make_pair(position_type(), structure_id_type());
};



#endif /* STRUCTUREFUNCTIONS_CPP */

