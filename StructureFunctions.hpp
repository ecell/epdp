#ifndef STRUCTFUNC_HPP
#define STRUCTFUNC_HPP

template <typename Ttraits_>
class PlanarSurface;
template <typename Ttraits_>
class CylindricalSurface;

template <typename Ttraits_>
inline std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>
get_pos(PlanarSurface<Ttraits_> const& origin_structure,
        PlanarSurface<Ttraits_> const& target_structure, typename Ttraits_::position_type const& old_pos);

template <typename Ttraits_>
inline std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>
get_pos(CylindricalSurface<Ttraits_> origin_structure,
        CylindricalSurface<Ttraits_> target_structure, typename Ttraits_::position_type const& old_pos);

#endif /* STRUCTUREFUNCTIONS_HPP */

