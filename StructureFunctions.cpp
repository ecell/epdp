#ifndef STRUCTFUNC_CPP
#define STRUCTFUNC_CPP

//#include "StructureFunctions.hpp" // TODO Why the hell does this not work, man???
#include "CuboidalRegion.hpp"
#include "SphericalSurface.hpp"
#include "CylindricalSurface.hpp"
#include "DiskSurface.hpp"
#include "PlanarSurface.hpp"


template <typename Ttraits_>
inline std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>
get_pos_sid_pair( CuboidalRegion<Ttraits_> const&         origin_structure,
                  CuboidalRegion<Ttraits_> const&         target_structure,
                  typename Ttraits_::position_type const& old_pos       )
{
    typedef typename Ttraits_::structure_id_type structure_id_type;
    typedef typename Ttraits_::position_type     position_type;

    return std::make_pair(position_type(), structure_id_type()); // TODO returns nothing meaningful for now
};

template <typename Ttraits_>
inline std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>
get_pos_sid_pair( SphericalSurface<Ttraits_> const&       origin_structure,
                  SphericalSurface<Ttraits_> const&       target_structure, 
                  typename Ttraits_::position_type const& old_pos       )
{
    typedef typename Ttraits_::structure_id_type structure_id_type;
    typedef typename Ttraits_::position_type     position_type;

    return std::make_pair(position_type(), structure_id_type()); // TODO returns nothing meaningful for now
};

template <typename Ttraits_>
inline std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>
get_pos_sid_pair( CylindricalSurface<Ttraits_> const&     origin_structure,
                  CylindricalSurface<Ttraits_> const&     target_structure, 
                  typename Ttraits_::position_type const& old_pos       )
{
    typedef typename Ttraits_::structure_id_type structure_id_type;
    typedef typename Ttraits_::position_type     position_type;

    return std::make_pair(position_type(), structure_id_type()); // TODO returns nothing meaningful for now
};

template <typename Ttraits_>
inline std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>
get_pos_sid_pair( CylindricalSurface<Ttraits_> const&     origin_structure,
                  DiskSurface<Ttraits_> const&            target_structure,
                  typename Ttraits_::position_type const& old_pos       )
{
    typedef typename Ttraits_::structure_id_type structure_id_type;
    typedef typename Ttraits_::position_type     position_type;

    return std::make_pair(position_type(), structure_id_type()); // TODO returns nothing meaningful for now
};

template <typename Ttraits_>
inline std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>
get_pos_sid_pair( PlanarSurface<Ttraits_> const&          origin_structure,
                  PlanarSurface<Ttraits_> const&          target_structure,
                  typename Ttraits_::position_type const& old_pos       )
{
    typedef typename Ttraits_::structure_id_type structure_id_type;
    typedef typename Ttraits_::position_type     position_type;

    const structure_id_type id(); // TODO what was this meant to be?
    return std::make_pair(target_structure.project_point(old_pos), target_structure.id());
};


#endif /* STRUCTUREFUNCTIONS_CPP */

