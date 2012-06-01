#ifndef STRUCTFUNC_CPP
#define STRUCTFUNC_CPP

#include "StructureFunctions.hpp"
#include "PlanarSurface.hpp"
#include "CylindricalSurface.hpp"
#include "DiskSurface.hpp"

template <typename Ttraits_>
inline std::pair<typename Ttraits::position_type, typename Ttraits::structure_id_type>
get_pos(PlanarSurface<Ttraits_> const& origin_structure,
        PlanarSurface<Ttraits_> const& target_structure, typename Ttraits_::position_type const& old_pos)
{
    typedef typename Ttraits_::structure_id_type structure_id_type;
    typedef typename Ttraits_::position_type     position_type;

    const structure_id_type id();
    return std::make_pair(target_structure.project_point(old_pos), target_structure.id());
}

template <typename Ttraits_>
inline std::pair<typename Ttraits::position_type, typename Ttraits::structure_id_type>
get_pos(CylindricalSurface<Ttraits_> origin_structure,
        CylindricalSurface<Ttraits_> target_structure, typename Ttraits_::position_type const& old_pos)
{
    typedef typename Ttraits_::structure_id_type structure_id_type;
    typedef typename Ttraits_::position_type     position_type;

    return std::make_pair(position_type(), structure_id_type())
}


#endif /* STRUCTUREFUNCTIONS_CPP */

