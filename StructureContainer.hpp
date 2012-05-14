#ifndef STRUCTURE_CONTAINER_HPP
#define STRUCTURE_CONTAINER_HPP

#include <boost/lexical_cast.hpp>
#include "exceptions.hpp"
#include "linear_algebra.hpp"
#include "ConnectivityContainer.hpp"
#include "CuboidalRegion.hpp"
#include "PlanarSurface.hpp"
#include "CylindricalSurface.hpp"

template <typename Tobj_, typename Tid_, typename Ttraits_>
class StructureContainer;

template<typename Ttraits_ >
std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>
apply_boundary (std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type> const& pos_structure_id,
                PlanarSurface<Ttraits_> const& planar_surface,
                StructureContainer<typename Ttraits_::structure_type, typename Ttraits_::structure_id_type, Ttraits_> const& sc);

template<typename Ttraits_ >
std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>
apply_boundary (std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type> const& pos_structure_id,
                CylindricalSurface<Ttraits_> const& cylindrical_surface,
                StructureContainer<typename Ttraits_::structure_type, typename Ttraits_::structure_id_type, Ttraits_> const& sc);


template <typename Tobj_, typename Tid_, typename Ttraits_>
class StructureContainer
{
public:
    typedef Tobj_               structure_type;
    typedef Tid_                structure_id_type;
    typedef Ttraits_            traits_type;
    typedef std::set<structure_id_type>                                             structure_id_set;
    typedef std::pair<const structure_id_type, boost::shared_ptr<structure_type> >  structure_id_pair;
    typedef std::map<structure_id_type, boost::shared_ptr<structure_type> >         structure_map;

protected:
    typedef select_second<typename structure_map::value_type>                       structures_second_selector_type;
    typedef boost::transform_iterator<structures_second_selector_type,
            typename structure_map::const_iterator>                                 structure_iterator;
    typedef typename std::map<structure_id_type, structure_id_set>                  per_structure_substructure_id_set;

    typedef CuboidalRegion<traits_type>                             cuboidalregion_type;
    typedef PlanarSurface<traits_type>                              planarsurface_type;
    typedef CylindricalSurface<traits_type>                         cylindricalsurface_type;

    typedef typename structure_type::position_type                      position_type;
    typedef typename structure_type::position_type                      vector_type;
    typedef std::pair<structure_id_type, vector_type>                   neighbor_id_vector_type;
    typedef ConnectivityContainer<structure_id_type, vector_type, 2>    cylinders_boundaryconditions_type;
    typedef ConnectivityContainer<structure_id_type, vector_type, 4>    planes_boundaryconditions_type;
    typedef ConnectivityContainer<structure_id_type, vector_type, 6>    boxes_boundaryconditions_type;

    typedef std::pair<structure_id_type, boost::shared_ptr<cuboidalregion_type> >       cuboidalreg_id_pair_type;
    typedef std::pair<structure_id_type, boost::shared_ptr<planarsurface_type> >        planarsurf_id_pair_type;
    typedef std::pair<structure_id_type, boost::shared_ptr<cylindricalsurface_type> >   cylindrsurf_id_pair_type;

public:
    typedef sized_iterator_range<structure_iterator>                                structures_range;
    typedef std::pair<position_type, structure_id_type>     position_structid_type;

public:
    virtual ~StructureContainer() {};

    // The constructor
    StructureContainer(structure_id_type default_structid) : default_structure_id_(default_structid)
    {
        structure_substructures_map_[default_structid] = structure_id_set();
    }

    virtual bool update_structure (cylindrsurf_id_pair_type const& structid_cylinder)  // can I make this into a template function?
    {
        if ( update_structure_base(structid_cylinder))
        {
            // add to Connectivity container for cylindrical surfaces
            cylindrical_structs_bc_.set_neighbor_info(structid_cylinder.first, 0, std::make_pair(structid_cylinder.first, multiply(structid_cylinder.second->shape().unit_z(), -1.0) ));
            cylindrical_structs_bc_.set_neighbor_info(structid_cylinder.first, 1, std::make_pair(structid_cylinder.first,          structid_cylinder.second->shape().unit_z() ));
            return true;
        }
        return false;
    }
    virtual bool update_structure (planarsurf_id_pair_type const& structid_plane)
    {
        // call regular update method (this also checks if the structure was already present)
        if ( update_structure_base(structid_plane))
        {
            // We now assume that the structure was not already in the boundary_conditions_thing
            // add to Connectivity container for planar surfaces and set reflective boundary conditions.
            planar_structs_bc_.set_neighbor_info(structid_plane.first, 0, std::make_pair(structid_plane.first, multiply(structid_plane.second->shape().unit_y(), -1.0) ));
            planar_structs_bc_.set_neighbor_info(structid_plane.first, 1, std::make_pair(structid_plane.first,          structid_plane.second->shape().unit_y() ));
            planar_structs_bc_.set_neighbor_info(structid_plane.first, 2, std::make_pair(structid_plane.first,          structid_plane.second->shape().unit_x() ));
            planar_structs_bc_.set_neighbor_info(structid_plane.first, 3, std::make_pair(structid_plane.first, multiply(structid_plane.second->shape().unit_x(), -1.0) ));
            return true;
        }
        // if the structure was already present, don't change anything
        return false;
    }
    virtual bool update_structure (cuboidalreg_id_pair_type const& structid_cube)
    {
        // add to Connectivity container for cuboidal regions
        // NOTE this information is never queried!!
        if ( update_structure_base(structid_cube))
        {
            cuboidal_structs_bc_.set_neighbor_info(structid_cube.first, 0, std::make_pair(structid_cube.first,          structid_cube.second->shape().unit_z()));
            cuboidal_structs_bc_.set_neighbor_info(structid_cube.first, 1, std::make_pair(structid_cube.first, multiply(structid_cube.second->shape().unit_z(), -1.0) ));
            cuboidal_structs_bc_.set_neighbor_info(structid_cube.first, 2, std::make_pair(structid_cube.first,          structid_cube.second->shape().unit_y()));
            cuboidal_structs_bc_.set_neighbor_info(structid_cube.first, 3, std::make_pair(structid_cube.first, multiply(structid_cube.second->shape().unit_y(), -1.0) ));
            cuboidal_structs_bc_.set_neighbor_info(structid_cube.first, 4, std::make_pair(structid_cube.first,          structid_cube.second->shape().unit_x()));
            cuboidal_structs_bc_.set_neighbor_info(structid_cube.first, 5, std::make_pair(structid_cube.first, multiply(structid_cube.second->shape().unit_x(), -1.0) ));
            return true;
        }
        return false;
    }
    virtual neighbor_id_vector_type get_neighbor_info(planarsurface_type const& structure, int n) const
    {
        planar_structs_bc_.get_neighbor_info(structure.id(), n);
    }
    virtual neighbor_id_vector_type get_neighbor_info(cylindricalsurface_type const& structure, int n) const
    {
        cylindrical_structs_bc_.get_neighbor_info(structure.id(), n);
    }
    template <typename Tstructure_>
    position_structid_type apply_boundary(Tstructure_ const& structure, position_structid_type const& pos_struct_id) const
    {
        return ::apply_boundary(pos_struct_id, structure, *this);
    }

    virtual structure_id_type get_def_structure_id() const
    {
        return default_structure_id_;
    }
    void set_def_structure(const boost::shared_ptr<cuboidalregion_type> cuboidalregion)
    {
        // check that the structure_type is the default structure_type
        if ( cuboidalregion->structure_id() != default_structure_id_)
        {
            throw illegal_state("Default structure should have itself as parent.");
        }
        cuboidalregion->set_id(default_structure_id_);
        update_structure(std::make_pair(default_structure_id_, cuboidalregion));
    }
    virtual bool has_structure(structure_id_type const& id) const
    {
        typename structure_map::const_iterator i(structure_map_.find(id));
        return i != structure_map_.end();
    }
    virtual boost::shared_ptr<structure_type> get_structure(structure_id_type const& id) const
    {
        typename structure_map::const_iterator i(structure_map_.find(id));
        if (structure_map_.end() == i)
        {
            throw not_found(std::string("Unknown structure (id=") + boost::lexical_cast<std::string>(id) + ")");
        }
        return (*i).second;
    }
    virtual structures_range get_structures_range() const
    {
        return structures_range(
        structure_iterator(structure_map_.begin(), structures_second_selector_type()),
        structure_iterator(structure_map_.end(),   structures_second_selector_type()),
        structure_map_.size());
    }
    virtual bool remove_structure(structure_id_type const& id)
    {
        // TODO
        //  -find StructureID in all ConnectivityContainers -> remove references
        //  -get all the substructures of structure
        //  -only remove if no substructures
        return false;
    }
    structure_id_set get_visible_structures(structure_id_type const& id) const
    {
        structure_id_set visible_structures;

        const boost::shared_ptr<const structure_type> structure(get_structure(id));
        if (structure->structure_id() == id)
        {
            visible_structures = structure_id_set();
        }
        else
        {
            visible_structures = get_visible_structures(structure->structure_id());
            visible_structures.erase(id);
        }
        const structure_id_set substructures (get_substructure_ids(id));
        visible_structures.insert(substructures.begin(), substructures.end());
        return visible_structures;
    }

    bool update_structure_base(structure_id_pair const& structid_pair)
    {
        typename structure_map::const_iterator i(structure_map_.find(structid_pair.first));
        if (i != structure_map_.end())
        // The item was already found in the map
        {
            if ((*i).second->structure_id() != structid_pair.second->structure_id())
            // If the structure had a different parent structure
            // Note that all the substructures of the structures are kept (they are aldo moved moved through the hierarchy)
            {
                structure_substructures_map_[(*i).second->structure_id()].erase((*i).first);
                structure_substructures_map_[structid_pair.second->structure_id()].insert(structid_pair.first);
            }
            structure_map_[(*i).first] = structid_pair.second;
            return false;
        }

        // The structure was not yet in the world.
        // create a new item in the structure mapping
        structure_map_[structid_pair.first] = structid_pair.second;
        // add the id the mapping 'super_structure_id -> set of substructures'
        structure_substructures_map_[structid_pair.second->structure_id()].insert(structid_pair.first);
        // create a new mapping from structure id -> set of substructures
        structure_substructures_map_[structid_pair.first] = structure_id_set();
        return true;
    }
private:
    // Get all the structure ids of the substructures
    structure_id_set get_substructure_ids(structure_id_type const& id) const
    {
        typename per_structure_substructure_id_set::const_iterator i(
            structure_substructures_map_.find(id));
        if (i == structure_substructures_map_.end())
        {
            throw not_found(std::string("Unknown structure (id=") + boost::lexical_cast<std::string>(id) + ")");
        }
        return (*i).second;
    }

protected:
    structure_map                       structure_map_;     // mapping: structure_id -> structure
    per_structure_substructure_id_set   structure_substructures_map_;
    structure_id_type                   default_structure_id_;

    // The next object describe the boundary conditions that are applied to the different types of structures.
    cylinders_boundaryconditions_type   cylindrical_structs_bc_;
    planes_boundaryconditions_type      planar_structs_bc_;
    boxes_boundaryconditions_type       cuboidal_structs_bc_;
};



//////// Inline functions applicable to the StructureContainer



//// functions for Planes
template<typename Ttraits_ >
inline std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>
apply_boundary (std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type> const& pos_structure_id,
                PlanarSurface<Ttraits_> const& planar_surface,
                StructureContainer<typename Ttraits_::structure_type, typename Ttraits_::structure_id_type, Ttraits_> const& sc)
// The template needs to be parameterized with the appropriate shape (which then parameterizes the type
// of the ConnectivityContainer that we use.
// We supply the structure container as an argument so that we can get the structures that we need.
// The connectivity container then needs to a part of the structure container such that we can also
// query the boundary conditions and connectivity.
{
    // useful typedefs
    typedef Ttraits_                    traits_type;
    typedef typename PlanarSurface<Ttraits_>::shape_type    plane_type;
//    typedef PlanarSurface<traits_type>  planar_surface_type;

    typedef typename plane_type::length_type            length_type;
    typedef typename plane_type::position_type          position_type;
    typedef typename traits_type::structure_id_type     structure_id_type;
    typedef std::pair<structure_id_type, position_type> neighbor_id_vector_type;

    // Note that we assume that the new position is in the plane (dot(pos, unit_z)==0)
//    const boost::shared_ptr<planar_surface_type> planar_surface(sc.get_structure(pos_structure_id.second));
    const plane_type origin_plane (planar_surface.shape());
    boost::array<length_type, 2> half_extends(origin_plane.half_extent());

    // TODO do cyclic transpose here first
    const position_type pos_vector( subtract(pos_structure_id.first, origin_plane.position()) );
    const length_type component_x ( dot_product(pos_vector, origin_plane.unit_x()) );
    const length_type component_y ( dot_product(pos_vector, origin_plane.unit_y()) );

    // declare the variables that will be written
    structure_id_type new_id(pos_structure_id.second);
    position_type neighbor_plane_par;
    position_type neighbor_plane_inl;

    if ( (-half_extends[0] < component_x) && ( component_x < half_extends[0]) &&
         (-half_extends[1] < component_y) && ( component_y < half_extends[1]) )
    {
        // we are still in the plane (did not pass any of the boundaries)
        // don't have to do anything
        return pos_structure_id;
    }
    else if ( (-half_extends[0] < component_x) && ( component_x < half_extends[0]) &&
              (half_extends[1] < component_y) )
    {
        // we are at the 'top' of the plane (side nr. 0)
        // get the unit vector pointing from the edge between the two planes to the center of the plane
        // TODO just pass origin_structure instead of plane and structure_id separately.
        const neighbor_id_vector_type neighbor_id_vector (sc.get_neighbor_info(planar_surface, 0));

        new_id = neighbor_id_vector.first;
        neighbor_plane_par = multiply(origin_plane.unit_x(), component_x);
        neighbor_plane_inl = add(multiply(origin_plane.unit_y(),     half_extends[1]),
                                 multiply(neighbor_id_vector.second, (component_y - half_extends[1]) ) );
    }
    else if ( (-half_extends[0] < component_x) && ( component_x < half_extends[0]) &&
              (component_y < -half_extends[1] ) )
    {
        // we are at the 'bottom' of the plane (side nr. 1)
        // get the unit vector pointing from the edge between the two planes to the center of the plane
        const neighbor_id_vector_type neighbor_id_vector (sc.get_neighbor_info(planar_surface, 1));

        new_id = neighbor_id_vector.first;
        neighbor_plane_par = multiply(origin_plane.unit_x(), component_x);
        neighbor_plane_inl = add(multiply(origin_plane.unit_y(),     -half_extends[1] ),
                                 multiply(neighbor_id_vector.second, (-component_y - half_extends[1]) ) );
    }
    else if ( (-half_extends[1] < component_y) && ( component_y < half_extends[1]) &&
              (component_x < -half_extends[0]) )
    {
        // we are at the 'left' of the plane (side nr. 2)
        // get the unit vector pointing from the edge between the two planes to the center of the plane
        const neighbor_id_vector_type neighbor_id_vector (sc.get_neighbor_info(planar_surface, 2));

        new_id = neighbor_id_vector.first;
        neighbor_plane_par = multiply(origin_plane.unit_y(), component_y);
        neighbor_plane_inl = add(multiply(origin_plane.unit_x(),     -half_extends[0]),
                                 multiply(neighbor_id_vector.second, (-component_x - half_extends[0]) ) );
    }
    else if ( (-half_extends[1] < component_y) && ( component_y < half_extends[1]) &&
              (half_extends[0] < component_x) )
    {
        // we are at the 'right' of the plane (side nr. 3)
        // this is a unit vector pointing from the edge between the two planes to the center of the plane
        const neighbor_id_vector_type neighbor_id_vector (sc.get_neighbor_info(planar_surface, 3));

        new_id = neighbor_id_vector.first;
        neighbor_plane_par = multiply(origin_plane.unit_y(), component_y);
        neighbor_plane_inl = add(multiply(origin_plane.unit_x(),     half_extends[0]),
                                 multiply(neighbor_id_vector.second, (component_x - half_extends[0]) ) );
    }
    else
    {
        // we are somewhere in the corners and don't know what to do
        // reject move?
        return pos_structure_id;
    }

    const position_type new_pos ( add(origin_plane.position(),
                                      add(neighbor_plane_par, neighbor_plane_inl)));
    return std::make_pair(new_pos, new_id);
}

///// Functions for Cylinders
template<typename Ttraits_ >
inline std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type>
apply_boundary (std::pair<typename Ttraits_::position_type, typename Ttraits_::structure_id_type> const& pos_structure_id,
                CylindricalSurface<Ttraits_> const& cylindrical_surface,
                StructureContainer<typename Ttraits_::structure_type, typename Ttraits_::structure_id_type, Ttraits_> const& sc)

// The template needs to be parameterized with the appropriate shape (which then parameterizes the type
// of the ConnectivityContainer that we use.
// We supply the structure container as an argument so that we can get the structures that we need.
// The connectivity container then needs to a part of the structure container such that we can also
// query the boundary conditions and connectivity.
{
    // useful typedefs
    typedef Ttraits_                        traits_type;
    typedef typename CylindricalSurface<Ttraits_>::shape_type    cylinder_type;
//    typedef CylindricalSurface<traits_type> cylindrical_surface_type;

    typedef typename cylinder_type::length_type         length_type;
    typedef typename cylinder_type::position_type       position_type;
    typedef typename traits_type::structure_id_type     structure_id_type;
    typedef std::pair<structure_id_type, position_type> neighbor_id_vector_type;

    // Note that we assume that the new position is in the cylinder (dot(pos, unit_r)==0)
//    const boost::shared_ptr<cylindrical_surface_type> cylindrical_surface(sc.get_structure(pos_structure_id.second));
    const cylinder_type origin_cylinder (cylindrical_surface.shape());
    const length_type half_length(origin_cylinder.half_length());

    // TODO do cyclic transpose here first
    const position_type pos_vector( subtract(pos_structure_id.first, origin_cylinder.position()) );
    const length_type component_z ( dot_product(pos_vector, origin_cylinder.unit_z()) );

    // declare the variables that will be written
    structure_id_type new_id(pos_structure_id.second);
    position_type neighbor_cylinder_inl;

    if ( half_length < component_z ) 
    {
        // we are at the 'right' side of the cylinder (side nr. 0, side in the direction of the normal vector)
        // get the unit vector pointing from the edge between the two planes to the center of the plane
        const neighbor_id_vector_type neighbor_id_vector (sc.get_neighbor_info(cylindrical_surface, 0));
//        const length_type half_length_new (sc.get_structure(new_id)->shape().half_length());

        new_id = neighbor_id_vector.first;
        neighbor_cylinder_inl = add( multiply(origin_cylinder.unit_z(),  half_length),
                                     multiply(neighbor_id_vector.second, (component_z - half_length) ) );
    }
    else if ( component_z < -half_length ) 
    {
        // we are at the 'right' side of the cylinder (side nr. 0, side in the direction of the normal vector)
        // get the unit vector pointing from the edge between the two planes to the center of the plane
        const neighbor_id_vector_type neighbor_id_vector (sc.get_neighbor_info(cylindrical_surface, 1));
//        const length_type half_length_new (sc.get_structure(new_id)->shape().half_length());

        new_id = neighbor_id_vector.first;
        neighbor_cylinder_inl = add( multiply(origin_cylinder.unit_z() , -half_length),
                                     multiply(neighbor_id_vector.second, (-component_z - half_length) ) );
    }
    else
    {
        // we are still in the cylinder (did not pass any of the boundaries)
        // don't have to do anything
        return pos_structure_id;
    }

    const position_type new_pos ( add(origin_cylinder.position(), neighbor_cylinder_inl) );
    return std::make_pair(new_pos, new_id);
}

#endif /* STRUCTURE_CONTAINER_HPP */

