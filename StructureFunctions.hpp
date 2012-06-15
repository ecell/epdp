#ifndef STRUCTFUNC_HPP
#define STRUCTFUNC_HPP


// Forward declaration of structures
template <typename Ttraits_>
class CuboidalRegion;

template <typename Ttraits_>
class SphericalSurface;

template <typename Ttraits_>
class CylindricalSurface;

template <typename Ttraits_>
class DiskSurface;

template <typename Ttraits_>
class PlanarSurface;


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
                  typename Ttraits_::position_type      const& old_pos          );

                  template <typename Ttraits_>
inline std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>
get_pos_sid_pair( CuboidalRegion<Ttraits_>              const& origin_structure,
                  SphericalSurface<Ttraits_>            const& target_structure,
                  typename Ttraits_::position_type      const& old_pos          );

template <typename Ttraits_>
inline std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>
get_pos_sid_pair( CuboidalRegion<Ttraits_>              const& origin_structure,
                  CylindricalSurface<Ttraits_>          const& target_structure,
                  typename Ttraits_::position_type      const& old_pos          );

template <typename Ttraits_>
inline std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>
get_pos_sid_pair( CuboidalRegion<Ttraits_>              const& origin_structure,
                  DiskSurface<Ttraits_>                 const& target_structure,
                  typename Ttraits_::position_type      const& old_pos          );

template <typename Ttraits_>
inline std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>
get_pos_sid_pair( CuboidalRegion<Ttraits_>              const& origin_structure,
                  PlanarSurface<Ttraits_>               const& target_structure,
                  typename Ttraits_::position_type      const& old_pos          );

/********************************/
/* Coming from SphericalSurface */
/********************************/
template <typename Ttraits_>
inline std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>
get_pos_sid_pair( SphericalSurface<Ttraits_>            const& origin_structure,
                  CuboidalRegion<Ttraits_>              const& target_structure,
                  typename Ttraits_::position_type      const& old_pos          );
                  
template <typename Ttraits_>
inline std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>
get_pos_sid_pair( SphericalSurface<Ttraits_>            const& origin_structure,
                  SphericalSurface<Ttraits_>            const& target_structure,
                  typename Ttraits_::position_type      const& old_pos          );
                  
template <typename Ttraits_>
inline std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>
get_pos_sid_pair( SphericalSurface<Ttraits_>            const& origin_structure,
                  CylindricalSurface<Ttraits_>          const& target_structure,
                  typename Ttraits_::position_type      const& old_pos          );
                  
template <typename Ttraits_>
inline std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>
get_pos_sid_pair( SphericalSurface<Ttraits_>            const& origin_structure,
                  DiskSurface<Ttraits_>                 const& target_structure,
                  typename Ttraits_::position_type      const& old_pos          );

template <typename Ttraits_>
inline std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>
get_pos_sid_pair( SphericalSurface<Ttraits_>            const& origin_structure,
                  PlanarSurface<Ttraits_>               const& target_structure,
                  typename Ttraits_::position_type      const& old_pos          );

/**********************************/
/* Coming from CylindricalSurface */
/**********************************/
template <typename Ttraits_>
inline std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>
get_pos_sid_pair( CylindricalSurface<Ttraits_>          const& origin_structure,
                  CuboidalRegion<Ttraits_>              const& target_structure,
                  typename Ttraits_::position_type      const& old_pos          );
                  
template <typename Ttraits_>
inline std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>
get_pos_sid_pair( CylindricalSurface<Ttraits_>          const& origin_structure,
                  SphericalSurface<Ttraits_>            const& target_structure,
                  typename Ttraits_::position_type      const& old_pos          );

template <typename Ttraits_>
inline std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>
get_pos_sid_pair( CylindricalSurface<Ttraits_>          const& origin_structure,
                  CylindricalSurface<Ttraits_>          const& target_structure,
                  typename Ttraits_::position_type      const& old_pos          );
                  
template <typename Ttraits_>
inline std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>
get_pos_sid_pair( CylindricalSurface<Ttraits_>          const& origin_structure,
                  DiskSurface<Ttraits_>                 const& target_structure, 
                  typename Ttraits_::position_type      const& old_pos          );
                  
template <typename Ttraits_>
inline std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>
get_pos_sid_pair( CylindricalSurface<Ttraits_>          const& origin_structure,
                  PlanarSurface<Ttraits_>               const& target_structure,
                  typename Ttraits_::position_type      const& old_pos          );

/***************************/
/* Coming from DiskSurface */
/***************************/
template <typename Ttraits_>
inline std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>
get_pos_sid_pair( DiskSurface<Ttraits_>                 const& origin_structure,
                  CuboidalRegion<Ttraits_>              const& target_structure, 
                  typename Ttraits_::position_type      const& old_pos          );
                  
template <typename Ttraits_>
inline std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>
get_pos_sid_pair( DiskSurface<Ttraits_>                 const& origin_structure,
                  SphericalSurface<Ttraits_>            const& target_structure, 
                  typename Ttraits_::position_type      const& old_pos          );
                  
template <typename Ttraits_>
inline std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>
get_pos_sid_pair( DiskSurface<Ttraits_>                 const& origin_structure,
                  CylindricalSurface<Ttraits_>          const& target_structure, 
                  typename Ttraits_::position_type      const& old_pos          );

template <typename Ttraits_>
inline std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>
get_pos_sid_pair( DiskSurface<Ttraits_>                 const& origin_structure,
                  DiskSurface<Ttraits_>                 const& target_structure, 
                  typename Ttraits_::position_type      const& old_pos          );
                  
template <typename Ttraits_>
inline std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>
get_pos_sid_pair( DiskSurface<Ttraits_>                 const& origin_structure,
                  PlanarSurface<Ttraits_>               const& target_structure, 
                  typename Ttraits_::position_type      const& old_pos          );
                  
/*****************************/
/* Coming from PlanarSurface */
/*****************************/
template <typename Ttraits_>
inline std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>
get_pos_sid_pair( PlanarSurface<Ttraits_>               const& origin_structure,
                  CuboidalRegion<Ttraits_>              const& target_structure, 
                  typename Ttraits_::position_type      const& old_pos          );
                  
template <typename Ttraits_>
inline std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>
get_pos_sid_pair( PlanarSurface<Ttraits_>               const& origin_structure,
                  SphericalSurface<Ttraits_>            const& target_structure, 
                  typename Ttraits_::position_type      const& old_pos          );
                  
template <typename Ttraits_>
inline std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>
get_pos_sid_pair( PlanarSurface<Ttraits_>               const& origin_structure,
                  CylindricalSurface<Ttraits_>          const& target_structure, 
                  typename Ttraits_::position_type      const& old_pos          );

template <typename Ttraits_>
inline std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>
get_pos_sid_pair( PlanarSurface<Ttraits_>               const& origin_structure,
                  DiskSurface<Ttraits_>                 const& target_structure, 
                  typename Ttraits_::position_type      const& old_pos          );
                  
template <typename Ttraits_>
inline std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>
get_pos_sid_pair( PlanarSurface<Ttraits_>               const& origin_structure,
                  PlanarSurface<Ttraits_>               const& target_structure, 
                  typename Ttraits_::position_type      const& old_pos          );

                  
#endif /* STRUCTUREFUNCTIONS_HPP */

