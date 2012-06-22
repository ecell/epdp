#ifndef STRUCTFUNC_CPP
#define STRUCTFUNC_CPP

//#include "StructureFunctions.hpp" // TODO Why the hell does this not work, man???
#include "CuboidalRegion.hpp"
#include "SphericalSurface.hpp"
#include "CylindricalSurface.hpp"
#include "DiskSurface.hpp"
#include "PlanarSurface.hpp"
#include "exceptions.hpp"


/******************************************************************************************/

// Define the structure functions and the actions performed in each case.
// There are two types: 1) yielding one new position and structure id
//                      2) yielding two new positions and structure ids
//                         (for single reactions with two products)

/******************************************************************************************/

/************************/
/*** ONE NEW POSITION ***/
/************************/

/******************************/
/* Coming from CuboidalRegion */
/******************************/
// CuboidalRegion -> CuboidalRegion
template <typename Ttraits_>
inline std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>
get_pos_sid_pair( CuboidalRegion<Ttraits_>              const& origin_structure,
                  CuboidalRegion<Ttraits_>              const& target_structure,
                  typename Ttraits_::position_type      const& old_pos          )
{
    typedef typename Ttraits_::structure_id_type        structure_id_type;
    typedef typename Ttraits_::position_type            position_type;
    typedef typename Ttraits_::length_type              length_type;
    
    // Currently only species change and decay are supported
    if(origin_structure.id == target_structure.id){
     
          return std::make_pair( old_pos, origin_structure.id );
    }
    else // structure transition not allowed
      
      throw illegal_propagation_attempt("Origin structure must be equal to target structure for this type of structure transition.");
};

// CuboidalRegion -> SphericalSurface
template <typename Ttraits_>
inline std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>
get_pos_sid_pair( CuboidalRegion<Ttraits_>              const& origin_structure,
                  SphericalSurface<Ttraits_>            const& target_structure,
                  typename Ttraits_::position_type      const& old_pos          )
{
    typedef typename Ttraits_::structure_id_type        structure_id_type;
    typedef typename Ttraits_::position_type            position_type;
    typedef typename Ttraits_::length_type              length_type;

    /*** COMBINATION NOT SUPPORTED ***/
    throw illegal_propagation_attempt("Structure transition between combination of origin structure and target structure not supported.");
    
    return std::make_pair(position_type(), structure_id_type());
};

// CuboidalRegion -> CylindricalSurface
template <typename Ttraits_>
inline std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>
get_pos_sid_pair( CuboidalRegion<Ttraits_>              const& origin_structure,
                  CylindricalSurface<Ttraits_>          const& target_structure,
                  typename Ttraits_::position_type      const& old_pos          )
{
    typedef typename Ttraits_::structure_id_type        structure_id_type;
    typedef typename Ttraits_::position_type            position_type;
    typedef typename Ttraits_::length_type              length_type;

    structure_id_type   new_id(    target_structure.id );
    position_type       new_pos(   target_structure->project_point(old_pos).first );
    length_type         proj_dist( target_structure->project_point(old_pos).second.second );
        // the distance of the projection of old_pos on target_structure to target_structure
    
    if(proj_dist < 0){ // if projection of old_pos is in structure
     
          return std::make_pair( new_pos, new_id );
    }
    else // structure transition not allowed
      
      throw illegal_propagation_attempt("Illegal original particle position for structure transition.");
    
};

// CuboidalRegion -> DiskSurface
template <typename Ttraits_>
inline std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>
get_pos_sid_pair( CuboidalRegion<Ttraits_>              const& origin_structure,
                  DiskSurface<Ttraits_>                 const& target_structure,
                  typename Ttraits_::position_type      const& old_pos          )
{
    typedef typename Ttraits_::structure_id_type        structure_id_type;
    typedef typename Ttraits_::position_type            position_type;
    typedef typename Ttraits_::length_type              length_type;

    /*** COMBINATION NOT SUPPORTED ***/
    throw illegal_propagation_attempt("Structure transition between combination of origin structure and target structure not supported.");
    
    return std::make_pair(position_type(), structure_id_type());
};

// CuboidalRegion -> PlanarSurface
template <typename Ttraits_>
inline std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>
get_pos_sid_pair( CuboidalRegion<Ttraits_>              const& origin_structure,
                  PlanarSurface<Ttraits_>               const& target_structure,
                  typename Ttraits_::position_type      const& old_pos          )
{
    typedef typename Ttraits_::structure_id_type        structure_id_type;
    typedef typename Ttraits_::position_type            position_type;
    typedef typename Ttraits_::length_type              length_type;

    structure_id_type   new_id(    target_structure.id );
    position_type       new_pos(   target_structure->project_point(old_pos).first );
    length_type         proj_dist( target_structure->project_point(old_pos).second.second );
        // the distance of the projection of old_pos on target_structure to target_structure
    
    if(proj_dist < 0){ // if projection of old_pos is in structure
     
          return std::make_pair( new_pos, new_id );
    }
    else // structure transition not allowed
      
      throw illegal_propagation_attempt("Illegal original particle position for structure transition.");
};

/********************************/
/* Coming from SphericalSurface */
/********************************/
// SphericalSurface -> CuboidalRegion
template <typename Ttraits_>
inline std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>
get_pos_sid_pair( SphericalSurface<Ttraits_>            const& origin_structure,
                  CuboidalRegion<Ttraits_>              const& target_structure,
                  typename Ttraits_::position_type      const& old_pos          )
{
    typedef typename Ttraits_::structure_id_type        structure_id_type;
    typedef typename Ttraits_::position_type            position_type;
    typedef typename Ttraits_::length_type              length_type;

    /*** COMBINATION NOT SUPPORTED ***/
    throw illegal_propagation_attempt("Structure transition between combination of origin structure and target structure not supported.");
    
    return std::make_pair(position_type(), structure_id_type());
};

// SphericalSurface -> SphericalSurface
template <typename Ttraits_>
inline std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>
get_pos_sid_pair( SphericalSurface<Ttraits_>            const& origin_structure,
                  SphericalSurface<Ttraits_>            const& target_structure,
                  typename Ttraits_::position_type      const& old_pos          )
{
    typedef typename Ttraits_::structure_id_type        structure_id_type;
    typedef typename Ttraits_::position_type            position_type;
    typedef typename Ttraits_::length_type              length_type;

    /*** COMBINATION NOT SUPPORTED ***/
    throw illegal_propagation_attempt("Structure transition between combination of origin structure and target structure not supported.");
    
    return std::make_pair(position_type(), structure_id_type());
};

// SphericalSurface -> CylindricalSurface
template <typename Ttraits_>
inline std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>
get_pos_sid_pair( SphericalSurface<Ttraits_>            const& origin_structure,
                  CylindricalSurface<Ttraits_>          const& target_structure,
                  typename Ttraits_::position_type      const& old_pos          )
{
    typedef typename Ttraits_::structure_id_type        structure_id_type;
    typedef typename Ttraits_::position_type            position_type;
    typedef typename Ttraits_::length_type              length_type;

    /*** COMBINATION NOT SUPPORTED ***/
    throw illegal_propagation_attempt("Structure transition between combination of origin structure and target structure not supported.");
    
    return std::make_pair(position_type(), structure_id_type());
};

// SphericalSurface -> DiskSurface
template <typename Ttraits_>
inline std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>
get_pos_sid_pair( SphericalSurface<Ttraits_>            const& origin_structure,
                  DiskSurface<Ttraits_>                 const& target_structure,
                  typename Ttraits_::position_type      const& old_pos          )
{
    typedef typename Ttraits_::structure_id_type        structure_id_type;
    typedef typename Ttraits_::position_type            position_type;
    typedef typename Ttraits_::length_type              length_type;

    /*** COMBINATION NOT SUPPORTED ***/
    throw illegal_propagation_attempt("Structure transition between combination of origin structure and target structure not supported.");
    
    return std::make_pair(position_type(), structure_id_type());
};

// SphericalSurface -> PlanarSurface
template <typename Ttraits_>
inline std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>
get_pos_sid_pair( SphericalSurface<Ttraits_>            const& origin_structure,
                  PlanarSurface<Ttraits_>               const& target_structure,
                  typename Ttraits_::position_type      const& old_pos          )
{
    typedef typename Ttraits_::structure_id_type        structure_id_type;
    typedef typename Ttraits_::position_type            position_type;
    typedef typename Ttraits_::length_type              length_type;

    /*** COMBINATION NOT SUPPORTED ***/
    throw illegal_propagation_attempt("Structure transition between combination of origin structure and target structure not supported.");
    
    return std::make_pair(position_type(), structure_id_type());
};

/**********************************/
/* Coming from CylindricalSurface */
/**********************************/
// CylindricalSurface -> CuboidalRegion
template <typename Ttraits_>
inline std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>
get_pos_sid_pair( CylindricalSurface<Ttraits_>          const& origin_structure,
                  CuboidalRegion<Ttraits_>              const& target_structure,
                  typename Ttraits_::position_type      const& old_pos          )
{
    typedef typename Ttraits_::structure_id_type        structure_id_type;
    typedef typename Ttraits_::position_type            position_type;
    typedef typename Ttraits_::length_type              length_type;

    /*** COMBINATION NOT SUPPORTED ***/
    throw illegal_propagation_attempt("Structure transition between combination of origin structure and target structure not supported.");
    
    return std::make_pair(position_type(), structure_id_type());
};

// CylindricalSurface -> CylindricalSurface
template <typename Ttraits_>
inline std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>
get_pos_sid_pair( CylindricalSurface<Ttraits_>          const& origin_structure,
                  CylindricalSurface<Ttraits_>          const& target_structure,
                  typename Ttraits_::position_type      const& old_pos          )
{
    typedef typename Ttraits_::structure_id_type        structure_id_type;
    typedef typename Ttraits_::position_type            position_type;
    typedef typename Ttraits_::length_type              length_type;

    // Currently only species change and decay are supported
    if(origin_structure.id == target_structure.id){
     
          return std::make_pair( old_pos, origin_structure.id );
    }
    else // structure transition not allowed
      
      throw illegal_propagation_attempt("Origin structure must be equal to target structure for this type of structure transition.");
};

// CylindricalSurface -> SphericalSurface
template <typename Ttraits_>
inline std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>
get_pos_sid_pair( CylindricalSurface<Ttraits_>          const& origin_structure,
                  SphericalSurface<Ttraits_>            const& target_structure,
                  typename Ttraits_::position_type      const& old_pos          )
{
    typedef typename Ttraits_::structure_id_type        structure_id_type;
    typedef typename Ttraits_::position_type            position_type;
    typedef typename Ttraits_::length_type              length_type;

    /*** COMBINATION NOT SUPPORTED ***/
    throw illegal_propagation_attempt("Structure transition between combination of origin structure and target structure not supported.");
    
    return std::make_pair(position_type(), structure_id_type());
};

// CylindricalSurface -> DiskSurface
template <typename Ttraits_>
inline std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>
get_pos_sid_pair( CylindricalSurface<Ttraits_>          const& origin_structure,
                  DiskSurface<Ttraits_>                 const& target_structure,
                  typename Ttraits_::position_type      const& old_pos          )
{
    typedef typename Ttraits_::structure_id_type        structure_id_type;
    typedef typename Ttraits_::position_type            position_type;
    typedef typename Ttraits_::length_type              length_type;

    structure_id_type   new_id(    target_structure.id );
    position_type       new_pos(   target_structure->project_point(old_pos).first );
    length_type         proj_dist( target_structure->project_point(old_pos).second.second );
        // the distance of the projection of old_pos on target_structure to target_structure
    
    if(proj_dist < 0){ // if projection of old_pos is in structure
    // TODO: Assert that DiskSurface is a substructure of CylindricalSurface ?
     
          return std::make_pair( new_pos, new_id );
    }
    else // structure transition not allowed
      
      throw illegal_propagation_attempt("Illegal original particle position for structure transition.");
};

// CylindricalSurface -> PlanarSurface
template <typename Ttraits_>
inline std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>
get_pos_sid_pair( CylindricalSurface<Ttraits_>          const& origin_structure,
                  PlanarSurface<Ttraits_>               const& target_structure,
                  typename Ttraits_::position_type      const& old_pos          )
{
    typedef typename Ttraits_::structure_id_type        structure_id_type;
    typedef typename Ttraits_::position_type            position_type;
    typedef typename Ttraits_::length_type              length_type;

    /*** COMBINATION NOT SUPPORTED ***/
    throw illegal_propagation_attempt("Structure transition between combination of origin structure and target structure not supported.");
    
    return std::make_pair(position_type(), structure_id_type());
};

/***************************/
/* Coming from DiskSurface */
/***************************/
// DiskSurface -> CuboidalRegion
template <typename Ttraits_>
inline std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>
get_pos_sid_pair( DiskSurface<Ttraits_>                 const& origin_structure,
                  CuboidalRegion<Ttraits_>              const& target_structure,
                  typename Ttraits_::position_type      const& old_pos          )
{
    typedef typename Ttraits_::structure_id_type        structure_id_type;
    typedef typename Ttraits_::position_type            position_type;
    typedef typename Ttraits_::length_type              length_type;

    /*** COMBINATION NOT SUPPORTED ***/
    throw illegal_propagation_attempt("Structure transition between combination of origin structure and target structure not supported.");
    
    return std::make_pair(position_type(), structure_id_type());
};

// DiskSurface -> SphericalSurface
template <typename Ttraits_>
inline std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>
get_pos_sid_pair( DiskSurface<Ttraits_>                 const& origin_structure,
                  SphericalSurface<Ttraits_>            const& target_structure,
                  typename Ttraits_::position_type      const& old_pos          )
{
    typedef typename Ttraits_::structure_id_type        structure_id_type;
    typedef typename Ttraits_::position_type            position_type;
    typedef typename Ttraits_::length_type              length_type;

    /*** COMBINATION NOT SUPPORTED ***/
    throw illegal_propagation_attempt("Structure transition between combination of origin structure and target structure not supported.");
    
    return std::make_pair(position_type(), structure_id_type());
};

// DiskSurface -> CylindricalSurface
template <typename Ttraits_>
inline std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>
get_pos_sid_pair( DiskSurface<Ttraits_>                 const& origin_structure,
                  CylindricalSurface<Ttraits_>          const& target_structure,
                  typename Ttraits_::position_type      const& old_pos          )
{
    typedef typename Ttraits_::structure_id_type        structure_id_type;
    typedef typename Ttraits_::position_type            position_type;
    typedef typename Ttraits_::length_type              length_type;

    /*** COMBINATION NOT SUPPORTED ***/
    throw illegal_propagation_attempt("Structure transition between combination of origin structure and target structure not supported.");
    
    return std::make_pair(position_type(), structure_id_type());
};

// DiskSurface -> DiskSurface
template <typename Ttraits_>
inline std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>
get_pos_sid_pair( DiskSurface<Ttraits_>                 const& origin_structure,
                  DiskSurface<Ttraits_>                 const& target_structure,
                  typename Ttraits_::position_type      const& old_pos          )
{
    typedef typename Ttraits_::structure_id_type        structure_id_type;
    typedef typename Ttraits_::position_type            position_type;
    typedef typename Ttraits_::length_type              length_type;

    // Currently only species change and decay are supported
    if(origin_structure.id == target_structure.id){
     
          return std::make_pair( old_pos, origin_structure.id );
    }
    else // structure transition not allowed
      
      throw illegal_propagation_attempt("Origin structure must be equal to target structure for this type of structure transition.");    
};

// DiskSurface -> PlanarSurface
template <typename Ttraits_>
inline std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>
get_pos_sid_pair( DiskSurface<Ttraits_>                 const& origin_structure,
                  PlanarSurface<Ttraits_>               const& target_structure,
                  typename Ttraits_::position_type      const& old_pos          )
{
    typedef typename Ttraits_::structure_id_type        structure_id_type;
    typedef typename Ttraits_::position_type            position_type;
    typedef typename Ttraits_::length_type              length_type;

    /*** COMBINATION NOT SUPPORTED ***/
    throw illegal_propagation_attempt("Structure transition between combination of origin structure and target structure not supported.");
    
    return std::make_pair(position_type(), structure_id_type());
};

/*****************************/
/* Coming from PlanarSurface */
/*****************************/
// PlanarSurface -> CuboidalRegion
template <typename Ttraits_>
inline std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>
get_pos_sid_pair( PlanarSurface<Ttraits_>               const& origin_structure,
                  CuboidalRegion<Ttraits_>              const& target_structure,
                  typename Ttraits_::position_type      const& old_pos          )
{
    typedef typename Ttraits_::structure_id_type        structure_id_type;
    typedef typename Ttraits_::position_type            position_type;
    typedef typename Ttraits_::length_type              length_type;

    /*** COMBINATION NOT SUPPORTED ***/
    throw illegal_propagation_attempt("Structure transition between combination of origin structure and target structure not supported.");
    
    return std::make_pair(position_type(), structure_id_type());
};

// PlanarSurface -> SphericalSurface
template <typename Ttraits_>
inline std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>
get_pos_sid_pair( PlanarSurface<Ttraits_>               const& origin_structure,
                  SphericalSurface<Ttraits_>            const& target_structure,
                  typename Ttraits_::position_type      const& old_pos          )
{
    typedef typename Ttraits_::structure_id_type        structure_id_type;
    typedef typename Ttraits_::position_type            position_type;
    typedef typename Ttraits_::length_type              length_type;

    /*** COMBINATION NOT SUPPORTED ***/
    throw illegal_propagation_attempt("Structure transition between combination of origin structure and target structure not supported.");
    
    return std::make_pair(position_type(), structure_id_type());
};

// PlanarSurface -> CylindricalSurface
template <typename Ttraits_>
inline std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>
get_pos_sid_pair( PlanarSurface<Ttraits_>               const& origin_structure,
                  CylindricalSurface<Ttraits_>          const& target_structure,
                  typename Ttraits_::position_type      const& old_pos          )
{
    typedef typename Ttraits_::structure_id_type        structure_id_type;
    typedef typename Ttraits_::position_type            position_type;
    typedef typename Ttraits_::length_type              length_type;

    /*** COMBINATION NOT SUPPORTED ***/
    throw illegal_propagation_attempt("Structure transition between combination of origin structure and target structure not supported.");
    
    return std::make_pair(position_type(), structure_id_type());
};

// PlanarSurface -> DiskSurface
template <typename Ttraits_>
inline std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>
get_pos_sid_pair( PlanarSurface<Ttraits_>               const& origin_structure,
                  DiskSurface<Ttraits_>                 const& target_structure,
                  typename Ttraits_::position_type      const& old_pos          )
{
    typedef typename Ttraits_::structure_id_type        structure_id_type;
    typedef typename Ttraits_::position_type            position_type;
    typedef typename Ttraits_::length_type              length_type;

    /*** COMBINATION NOT SUPPORTED ***/
    throw illegal_propagation_attempt("Structure transition between combination of origin structure and target structure not supported.");
    
    return std::make_pair(position_type(), structure_id_type());
};

// PlanarSurface -> PlanarSurface
template <typename Ttraits_>
inline std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>
get_pos_sid_pair( PlanarSurface<Ttraits_>               const& origin_structure,
                  PlanarSurface<Ttraits_>               const& target_structure,
                  typename Ttraits_::position_type      const& old_pos          )
{
    typedef typename Ttraits_::structure_id_type        structure_id_type;
    typedef typename Ttraits_::position_type            position_type;
    typedef typename Ttraits_::length_type              length_type;

    // Currently only species change and decay are supported
    // TODO Implement transitions from one plane to another
    if(origin_structure.id == target_structure.id){
     
          return std::make_pair( old_pos, origin_structure.id );
    }
    else // structure transition not allowed
      
      throw illegal_propagation_attempt("Origin structure must be equal to target structure for this type of structure transition.");        
};


/******************************************************************************************/


/*************************/
/*** TWO NEW POSITIONS ***/
/*************************/

/******************************/
/* Coming from CuboidalRegion */
/******************************/
// CuboidalRegion -> CuboidalRegion
template <typename Ttraits_>
inline std::pair< std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>, 
                  std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type> >
get_pos_sid_pair_pair( CuboidalRegion<Ttraits_>              const& origin_structure,
                       CuboidalRegion<Ttraits_>              const& target_structure,
                       typename Ttraits_::position_type      const& old_pos          )
{
    typedef typename Ttraits_::structure_id_type        structure_id_type;
    typedef typename Ttraits_::position_type            position_type;
    typedef typename Ttraits_::length_type              length_type;

    /*** COMBINATION NOT SUPPORTED ***/
    throw illegal_propagation_attempt("Structure transition between combination of origin structure and target structure not supported.");
    
    return std::make_pair(      std::make_pair(position_type(), structure_id_type()),
                                std::make_pair(position_type(), structure_id_type())    );
};

// CuboidalRegion -> SphericalSurface
template <typename Ttraits_>
inline std::pair< std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>, 
                  std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type> >
get_pos_sid_pair_pair( CuboidalRegion<Ttraits_>              const& origin_structure,
                       SphericalSurface<Ttraits_>            const& target_structure,
                       typename Ttraits_::position_type      const& old_pos          )
{
    typedef typename Ttraits_::structure_id_type        structure_id_type;
    typedef typename Ttraits_::position_type            position_type;
    typedef typename Ttraits_::length_type              length_type;

    /*** COMBINATION NOT SUPPORTED ***/
    throw illegal_propagation_attempt("Structure transition between combination of origin structure and target structure not supported.");
    
    return std::make_pair(      std::make_pair(position_type(), structure_id_type()),
                                std::make_pair(position_type(), structure_id_type())    );
};

// CuboidalRegion -> CylindricalSurface
template <typename Ttraits_>
inline std::pair< std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>, 
                  std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type> >
get_pos_sid_pair_pair( CuboidalRegion<Ttraits_>              const& origin_structure,
                       CylindricalSurface<Ttraits_>          const& target_structure,
                       typename Ttraits_::position_type      const& old_pos          )
{
    typedef typename Ttraits_::structure_id_type        structure_id_type;
    typedef typename Ttraits_::position_type            position_type;
    typedef typename Ttraits_::length_type              length_type;

    /*** COMBINATION NOT SUPPORTED ***/
    throw illegal_propagation_attempt("Structure transition between combination of origin structure and target structure not supported.");
    
    return std::make_pair(      std::make_pair(position_type(), structure_id_type()),
                                std::make_pair(position_type(), structure_id_type())    );
    
};

// CuboidalRegion -> DiskSurface
template <typename Ttraits_>
inline std::pair< std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>, 
                  std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type> >
get_pos_sid_pair_pair( CuboidalRegion<Ttraits_>              const& origin_structure,
                       DiskSurface<Ttraits_>                 const& target_structure,
                       typename Ttraits_::position_type      const& old_pos          )
{
    typedef typename Ttraits_::structure_id_type        structure_id_type;
    typedef typename Ttraits_::position_type            position_type;
    typedef typename Ttraits_::length_type              length_type;

    /*** COMBINATION NOT SUPPORTED ***/
    throw illegal_propagation_attempt("Structure transition between combination of origin structure and target structure not supported.");
    
    return std::make_pair(      std::make_pair(position_type(), structure_id_type()),
                                std::make_pair(position_type(), structure_id_type())    );
};

// CuboidalRegion -> PlanarSurface
template <typename Ttraits_>
inline std::pair< std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>, 
                  std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type> >
get_pos_sid_pair_pair( CuboidalRegion<Ttraits_>              const& origin_structure,
                       PlanarSurface<Ttraits_>               const& target_structure,
                       typename Ttraits_::position_type      const& old_pos          )
{
    typedef typename Ttraits_::structure_id_type        structure_id_type;
    typedef typename Ttraits_::position_type            position_type;
    typedef typename Ttraits_::length_type              length_type;

    /*** COMBINATION NOT SUPPORTED ***/
    throw illegal_propagation_attempt("Structure transition between combination of origin structure and target structure not supported.");
    
    return std::make_pair(      std::make_pair(position_type(), structure_id_type()),
                                std::make_pair(position_type(), structure_id_type())    );
};

/********************************/
/* Coming from SphericalSurface */
/********************************/
// SphericalSurface -> CuboidalRegion
template <typename Ttraits_>
inline std::pair< std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>, 
                  std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type> >
get_pos_sid_pair_pair( SphericalSurface<Ttraits_>            const& origin_structure,
                       CuboidalRegion<Ttraits_>              const& target_structure,
                       typename Ttraits_::position_type      const& old_pos          )
{
    typedef typename Ttraits_::structure_id_type        structure_id_type;
    typedef typename Ttraits_::position_type            position_type;
    typedef typename Ttraits_::length_type              length_type;

    /*** COMBINATION NOT SUPPORTED ***/
    throw illegal_propagation_attempt("Structure transition between combination of origin structure and target structure not supported.");
    
    return std::make_pair(      std::make_pair(position_type(), structure_id_type()),
                                std::make_pair(position_type(), structure_id_type())    );
};

// SphericalSurface -> SphericalSurface
template <typename Ttraits_>
inline std::pair< std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>, 
                  std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type> >
get_pos_sid_pair_pair( SphericalSurface<Ttraits_>            const& origin_structure,
                       SphericalSurface<Ttraits_>            const& target_structure,
                       typename Ttraits_::position_type      const& old_pos          )
{
    typedef typename Ttraits_::structure_id_type        structure_id_type;
    typedef typename Ttraits_::position_type            position_type;
    typedef typename Ttraits_::length_type              length_type;

    /*** COMBINATION NOT SUPPORTED ***/
    throw illegal_propagation_attempt("Structure transition between combination of origin structure and target structure not supported.");
    
    return std::make_pair(      std::make_pair(position_type(), structure_id_type()),
                                std::make_pair(position_type(), structure_id_type())    );
};

// SphericalSurface -> CylindricalSurface
template <typename Ttraits_>
inline std::pair< std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>, 
                  std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type> >
get_pos_sid_pair_pair( SphericalSurface<Ttraits_>            const& origin_structure,
                       CylindricalSurface<Ttraits_>          const& target_structure,
                       typename Ttraits_::position_type      const& old_pos          )
{
    typedef typename Ttraits_::structure_id_type        structure_id_type;
    typedef typename Ttraits_::position_type            position_type;
    typedef typename Ttraits_::length_type              length_type;

    /*** COMBINATION NOT SUPPORTED ***/
    throw illegal_propagation_attempt("Structure transition between combination of origin structure and target structure not supported.");
    
    return std::make_pair(      std::make_pair(position_type(), structure_id_type()),
                                std::make_pair(position_type(), structure_id_type())    );
};

// SphericalSurface -> DiskSurface
template <typename Ttraits_>
inline std::pair< std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>, 
                  std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type> >
get_pos_sid_pair_pair( SphericalSurface<Ttraits_>            const& origin_structure,
                       DiskSurface<Ttraits_>                 const& target_structure,
                       typename Ttraits_::position_type      const& old_pos          )
{
    typedef typename Ttraits_::structure_id_type        structure_id_type;
    typedef typename Ttraits_::position_type            position_type;
    typedef typename Ttraits_::length_type              length_type;

    /*** COMBINATION NOT SUPPORTED ***/
    throw illegal_propagation_attempt("Structure transition between combination of origin structure and target structure not supported.");
    
    return std::make_pair(      std::make_pair(position_type(), structure_id_type()),
                                std::make_pair(position_type(), structure_id_type())    );
};

// SphericalSurface -> PlanarSurface
template <typename Ttraits_>
inline std::pair< std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>, 
                  std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type> >
get_pos_sid_pair_pair( SphericalSurface<Ttraits_>            const& origin_structure,
                       PlanarSurface<Ttraits_>               const& target_structure,
                       typename Ttraits_::position_type      const& old_pos          )
{
    typedef typename Ttraits_::structure_id_type        structure_id_type;
    typedef typename Ttraits_::position_type            position_type;
    typedef typename Ttraits_::length_type              length_type;

    /*** COMBINATION NOT SUPPORTED ***/
    throw illegal_propagation_attempt("Structure transition between combination of origin structure and target structure not supported.");
    
    return std::make_pair(      std::make_pair(position_type(), structure_id_type()),
                                std::make_pair(position_type(), structure_id_type())    );
};

/**********************************/
/* Coming from CylindricalSurface */
/**********************************/
// CylindricalSurface -> CuboidalRegion
template <typename Ttraits_>
inline std::pair< std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>, 
                  std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type> >
get_pos_sid_pair_pair( CylindricalSurface<Ttraits_>          const& origin_structure,
                       CuboidalRegion<Ttraits_>              const& target_structure,
                       typename Ttraits_::position_type      const& old_pos          )
{
    typedef typename Ttraits_::structure_id_type        structure_id_type;
    typedef typename Ttraits_::position_type            position_type;
    typedef typename Ttraits_::length_type              length_type;

    /*** COMBINATION NOT SUPPORTED ***/
    throw illegal_propagation_attempt("Structure transition between combination of origin structure and target structure not supported.");
    
    return std::make_pair(      std::make_pair(position_type(), structure_id_type()),
                                std::make_pair(position_type(), structure_id_type())    );
};

// CylindricalSurface -> CylindricalSurface
template <typename Ttraits_>
inline std::pair< std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>, 
                  std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type> >
get_pos_sid_pair_pair( CylindricalSurface<Ttraits_>          const& origin_structure,
                       CylindricalSurface<Ttraits_>          const& target_structure,
                       typename Ttraits_::position_type      const& old_pos          )
{
    typedef typename Ttraits_::structure_id_type        structure_id_type;
    typedef typename Ttraits_::position_type            position_type;
    typedef typename Ttraits_::length_type              length_type;

    /*** COMBINATION NOT SUPPORTED ***/
    throw illegal_propagation_attempt("Structure transition between combination of origin structure and target structure not supported.");
    
    return std::make_pair(      std::make_pair(position_type(), structure_id_type()),
                                std::make_pair(position_type(), structure_id_type())    );
};

// CylindricalSurface -> SphericalSurface
template <typename Ttraits_>
inline std::pair< std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>, 
                  std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type> >
get_pos_sid_pair_pair( CylindricalSurface<Ttraits_>          const& origin_structure,
                       SphericalSurface<Ttraits_>            const& target_structure,
                       typename Ttraits_::position_type      const& old_pos          )
{
    typedef typename Ttraits_::structure_id_type        structure_id_type;
    typedef typename Ttraits_::position_type            position_type;
    typedef typename Ttraits_::length_type              length_type;

    /*** COMBINATION NOT SUPPORTED ***/
    throw illegal_propagation_attempt("Structure transition between combination of origin structure and target structure not supported.");
    
    return std::make_pair(      std::make_pair(position_type(), structure_id_type()),
                                std::make_pair(position_type(), structure_id_type())    );
};

// CylindricalSurface -> DiskSurface
template <typename Ttraits_>
inline std::pair< std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>, 
                  std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type> >
get_pos_sid_pair_pair( CylindricalSurface<Ttraits_>          const& origin_structure,
                       DiskSurface<Ttraits_>                 const& target_structure,
                       typename Ttraits_::position_type      const& old_pos          )
{
    typedef typename Ttraits_::structure_id_type        structure_id_type;
    typedef typename Ttraits_::position_type            position_type;
    typedef typename Ttraits_::length_type              length_type;

    /*** COMBINATION NOT SUPPORTED ***/
    throw illegal_propagation_attempt("Structure transition between combination of origin structure and target structure not supported.");
    
    return std::make_pair(      std::make_pair(position_type(), structure_id_type()),
                                std::make_pair(position_type(), structure_id_type())    );
};

// CylindricalSurface -> PlanarSurface
template <typename Ttraits_>
inline std::pair< std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>, 
                  std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type> >
get_pos_sid_pair_pair( CylindricalSurface<Ttraits_>          const& origin_structure,
                       PlanarSurface<Ttraits_>               const& target_structure,
                       typename Ttraits_::position_type      const& old_pos          )
{
    typedef typename Ttraits_::structure_id_type        structure_id_type;
    typedef typename Ttraits_::position_type            position_type;
    typedef typename Ttraits_::length_type              length_type;

    /*** COMBINATION NOT SUPPORTED ***/
    throw illegal_propagation_attempt("Structure transition between combination of origin structure and target structure not supported.");
    
    return std::make_pair(      std::make_pair(position_type(), structure_id_type()),
                                std::make_pair(position_type(), structure_id_type())    );
};

/***************************/
/* Coming from DiskSurface */
/***************************/
// DiskSurface -> CuboidalRegion
template <typename Ttraits_>
inline std::pair< std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>, 
                  std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type> >
get_pos_sid_pair_pair( DiskSurface<Ttraits_>                 const& origin_structure,
                       CuboidalRegion<Ttraits_>              const& target_structure,
                       typename Ttraits_::position_type      const& old_pos          )
{
    typedef typename Ttraits_::structure_id_type        structure_id_type;
    typedef typename Ttraits_::position_type            position_type;
    typedef typename Ttraits_::length_type              length_type;

    /*** COMBINATION NOT SUPPORTED ***/
    throw illegal_propagation_attempt("Structure transition between combination of origin structure and target structure not supported.");
    
    return std::make_pair(      std::make_pair(position_type(), structure_id_type()),
                                std::make_pair(position_type(), structure_id_type())    );
};

// DiskSurface -> SphericalSurface
template <typename Ttraits_>
inline std::pair< std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>, 
                  std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type> >
get_pos_sid_pair_pair( DiskSurface<Ttraits_>                 const& origin_structure,
                       SphericalSurface<Ttraits_>            const& target_structure,
                       typename Ttraits_::position_type      const& old_pos          )
{
    typedef typename Ttraits_::structure_id_type        structure_id_type;
    typedef typename Ttraits_::position_type            position_type;
    typedef typename Ttraits_::length_type              length_type;

    /*** COMBINATION NOT SUPPORTED ***/
    throw illegal_propagation_attempt("Structure transition between combination of origin structure and target structure not supported.");
    
    return std::make_pair(      std::make_pair(position_type(), structure_id_type()),
                                std::make_pair(position_type(), structure_id_type())    );
};

// DiskSurface -> CylindricalSurface
template <typename Ttraits_>
inline std::pair< std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>, 
                  std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type> >
get_pos_sid_pair_pair( DiskSurface<Ttraits_>                 const& origin_structure,
                       CylindricalSurface<Ttraits_>          const& target_structure,
                       typename Ttraits_::position_type      const& old_pos          )
{
    typedef typename Ttraits_::structure_id_type        structure_id_type;
    typedef typename Ttraits_::position_type            position_type;
    typedef typename Ttraits_::length_type              length_type;

    /*** COMBINATION NOT SUPPORTED ***/
    throw illegal_propagation_attempt("Structure transition between combination of origin structure and target structure not supported.");
    
    return std::make_pair(      std::make_pair(position_type(), structure_id_type()),
                                std::make_pair(position_type(), structure_id_type())    );
};

// DiskSurface -> DiskSurface
template <typename Ttraits_>
inline std::pair< std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>, 
                  std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type> >
get_pos_sid_pair_pair( DiskSurface<Ttraits_>                 const& origin_structure,
                       DiskSurface<Ttraits_>                 const& target_structure,
                       typename Ttraits_::position_type      const& old_pos          )
{
    typedef typename Ttraits_::structure_id_type        structure_id_type;
    typedef typename Ttraits_::position_type            position_type;
    typedef typename Ttraits_::length_type              length_type;

    /*** COMBINATION NOT SUPPORTED ***/
    throw illegal_propagation_attempt("Structure transition between combination of origin structure and target structure not supported.");
    
    return std::make_pair(      std::make_pair(position_type(), structure_id_type()),
                                std::make_pair(position_type(), structure_id_type())    );
};

// DiskSurface -> PlanarSurface
template <typename Ttraits_>
inline std::pair< std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>, 
                  std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type> >
get_pos_sid_pair_pair( DiskSurface<Ttraits_>                 const& origin_structure,
                       PlanarSurface<Ttraits_>               const& target_structure,
                       typename Ttraits_::position_type      const& old_pos          )
{
    typedef typename Ttraits_::structure_id_type        structure_id_type;
    typedef typename Ttraits_::position_type            position_type;
    typedef typename Ttraits_::length_type              length_type;

    /*** COMBINATION NOT SUPPORTED ***/
    throw illegal_propagation_attempt("Structure transition between combination of origin structure and target structure not supported.");
    
    return std::make_pair(      std::make_pair(position_type(), structure_id_type()),
                                std::make_pair(position_type(), structure_id_type())    );
};

/*****************************/
/* Coming from PlanarSurface */
/*****************************/
// PlanarSurface -> CuboidalRegion
template <typename Ttraits_>
inline std::pair< std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>, 
                  std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type> >
get_pos_sid_pair_pair( PlanarSurface<Ttraits_>               const& origin_structure,
                       CuboidalRegion<Ttraits_>              const& target_structure,
                       typename Ttraits_::position_type      const& old_pos          )
{
    typedef typename Ttraits_::structure_id_type        structure_id_type;
    typedef typename Ttraits_::position_type            position_type;
    typedef typename Ttraits_::length_type              length_type;

    /*** COMBINATION NOT SUPPORTED ***/
    throw illegal_propagation_attempt("Structure transition between combination of origin structure and target structure not supported.");
    
    return std::make_pair(      std::make_pair(position_type(), structure_id_type()),
                                std::make_pair(position_type(), structure_id_type())    );
};

// PlanarSurface -> SphericalSurface
template <typename Ttraits_>
inline std::pair< std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>, 
                  std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type> >
get_pos_sid_pair_pair( PlanarSurface<Ttraits_>               const& origin_structure,
                       SphericalSurface<Ttraits_>            const& target_structure,
                       typename Ttraits_::position_type      const& old_pos          )
{
    typedef typename Ttraits_::structure_id_type        structure_id_type;
    typedef typename Ttraits_::position_type            position_type;
    typedef typename Ttraits_::length_type              length_type;

    /*** COMBINATION NOT SUPPORTED ***/
    throw illegal_propagation_attempt("Structure transition between combination of origin structure and target structure not supported.");
    
    return std::make_pair(      std::make_pair(position_type(), structure_id_type()),
                                std::make_pair(position_type(), structure_id_type())    );
};

// PlanarSurface -> CylindricalSurface
template <typename Ttraits_>
inline std::pair< std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>, 
                  std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type> >
get_pos_sid_pair_pair( PlanarSurface<Ttraits_>               const& origin_structure,
                       CylindricalSurface<Ttraits_>          const& target_structure,
                       typename Ttraits_::position_type      const& old_pos          )
{
    typedef typename Ttraits_::structure_id_type        structure_id_type;
    typedef typename Ttraits_::position_type            position_type;
    typedef typename Ttraits_::length_type              length_type;

    /*** COMBINATION NOT SUPPORTED ***/
    throw illegal_propagation_attempt("Structure transition between combination of origin structure and target structure not supported.");
    
    return std::make_pair(      std::make_pair(position_type(), structure_id_type()),
                                std::make_pair(position_type(), structure_id_type())    );
};

// PlanarSurface -> DiskSurface
template <typename Ttraits_>
inline std::pair< std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>, 
                  std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type> >
get_pos_sid_pair_pair( PlanarSurface<Ttraits_>               const& origin_structure,
                       DiskSurface<Ttraits_>                 const& target_structure,
                       typename Ttraits_::position_type      const& old_pos          )
{
    typedef typename Ttraits_::structure_id_type        structure_id_type;
    typedef typename Ttraits_::position_type            position_type;
    typedef typename Ttraits_::length_type              length_type;

    /*** COMBINATION NOT SUPPORTED ***/
    throw illegal_propagation_attempt("Structure transition between combination of origin structure and target structure not supported.");
    
    return std::make_pair(      std::make_pair(position_type(), structure_id_type()),
                                std::make_pair(position_type(), structure_id_type())    );
};

// PlanarSurface -> PlanarSurface
template <typename Ttraits_>
inline std::pair< std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>, 
                  std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type> >
get_pos_sid_pair_pair( PlanarSurface<Ttraits_>               const& origin_structure,
                       PlanarSurface<Ttraits_>               const& target_structure,
                       typename Ttraits_::position_type      const& old_pos          )
{
    typedef typename Ttraits_::structure_id_type        structure_id_type;
    typedef typename Ttraits_::position_type            position_type;
    typedef typename Ttraits_::length_type              length_type;

    /*** COMBINATION NOT SUPPORTED ***/
    throw illegal_propagation_attempt("Structure transition between combination of origin structure and target structure not supported.");
    
    return std::make_pair(      std::make_pair(position_type(), structure_id_type()),
                                std::make_pair(position_type(), structure_id_type())    );
};

/******************************************************************************************************/


#endif /* STRUCTUREFUNCTIONS_CPP */

