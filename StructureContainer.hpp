#ifndef STRUCTURE_CONTAINER_HPP
#define STRUCTURE_CONTAINER_HPP

#include <boost/lexical_cast.hpp>
#include "exceptions.hpp"
#include "linear_algebra.hpp"

/*
template<typename Ttraits, Tobj_, std::size_t max_neighbors>
class ConnectivityContainer
// implements the way structures are connected.
// It allows you to check what structure you enter into when you leave the current one.
{
public:
    typedef typename Ttraits        traits_type;
    typedef typename Tobj_          obj_type;
    typedef typename std::size_t    size_type;
    typedef typename traits_type::structure_id_type structure_id_type;
    typedef typename traits_type::position_type     position_type;

protected:
    // the structureID of the neighbor and the normal vector pointing from the
    // intersection in the direction of the center of the neighbor
    typedef std::pair<structure_id_type, position_type>                     structure_id_vector_pair_type;
    typedef boost::array<structure_id_vector_pair_type, max_neighbors>      structure_id_vector_array_type;
    // the mapping from the structure and all its neighbors and how they are connected.
    typedef std::map<structure_id_type, structure_id_vector_array_type>     structid_neighbor_pair_map;
    typedef select_second<typename structid_neighbor_pair_map::value_type>  structid_neighbor_second_selector_type;
    typedef boost::transform_iterator<structid_neighbor_second_selector_type,
            typename structid_neighbor_pair_map::const_iterator>            structid_neighbor_iterator;

public:
    static const size_type max = max_neighbors;

public:
    virtual bool add_neighbor_info(obj_type const& obj) = 0;

    virtual structure_id_vector_pair_type get_neighbor_info(structure_id_type const& id, size_type const n) const
    {
        // check that the index is not too large
        if ( n>= max )
            throw not_found(std::string("Index out of range for neighbor (n= ") + boost::lexical_cast<std::string>(n) + ")");

        typename structid_neighbor_pair_map::const_iterator i(structid_neighbor_pair_map_.find(id));
        if (structid_neighbor_pair_map_.end() == i)
        {
            throw not_found(std::string("Neighbor structure not found(id=") + boost::lexical_cast<std::string>(id) + ")");
        }
        return (*i).second[n];
    }

    virtual ConnectivityContainer()
    // do we start with an empty container or do we construct it as a static object with n structures already connected.

private:
    // mapping StructureID -> array<(StructureID, vector3D), max_neighbors>
    structid_neighbor_pair_map  neighbor_mapping_;
};

*/




template <typename Tobj_, typename Tid_>
class StructureContainer
{
public:
    typedef Tobj_              structure_type;
    typedef Tid_               structure_id_type;
    typedef std::set<structure_id_type>                                             structure_id_set;
    typedef std::pair<const structure_id_type, boost::shared_ptr<structure_type> >  structure_id_pair;
    typedef std::map<structure_id_type, boost::shared_ptr<structure_type> >         structure_map;

protected:
    typedef select_second<typename structure_map::value_type>                       structures_second_selector_type;
    typedef boost::transform_iterator<structures_second_selector_type,
            typename structure_map::const_iterator>                                 structure_iterator;
    typedef typename std::map<structure_id_type, structure_id_set>                  per_structure_substructure_id_set;

public:
    typedef sized_iterator_range<structure_iterator>                                structures_range;


public:
    virtual ~StructureContainer() {};

    StructureContainer(structure_id_type default_structid) : default_structure_id_(default_structid)
    {
        structure_substructures_map_[default_structid] = structure_id_set();
    }

/*    add_structure (CylindricalSurface)  // can I make this into a template function?
    {
        // add to Connectivity container for cylindrical surfaces
    }
    add_structure (PlanarSurface)
    {
        // add to Connectivity container for planar surfaces
    }
    add_structure (CuboidalRegion)
    {
        // add to Connectivity container for planar surfaces
    }
    remove_structure (StructureID)
    {
        // find StructureID in all ConnectivityContainers -> remove references
    }

*/
    virtual structure_id_type get_def_structure_id() const
    {
        return default_structure_id_;
    }
    void set_def_structure(const boost::shared_ptr<structure_type> structure)
    {
        // check that the structure_type is the default structure_type
        if ( structure->structure_id() != default_structure_id_)
        {
            throw illegal_state("Default structure doesn't have right properties");
        }
        structure->set_id(default_structure_id_);
        update_structure(std::make_pair(default_structure_id_, structure));
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
    virtual bool update_structure(structure_id_pair const& structid_pair)
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
    virtual bool remove_structure(structure_id_type const& id)
    {
        // TODO
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

//    CylindricalSurfaceConnectivityContainer     // This contains a StructureID -> vector <(StructureID, index), 2> map
//    PlanarSurfaceConnectivityContainer          // This contains a StructureID -> vector <(StructureID, index), 4> map
//    CuboidalReginConnectivityContainer          // This contains a StructureID -> vector <(StructureID, index), 6> map

};
/*


//////// Inline functions applicable to the StructureContainer


//// functions for Planes

template<typename ConnectivityContainer<> >
inline typename ConnectivityContainer::traits_type::position_structid_pair_type
apply_boundary (StructureContainer const& sc, position_structid_pair_type const& pos_structure_id) const
// The template needs to be parameterized with the appropriate shape (which then parameterizes the type
// of the ConnectivityContainer that we use.
// We supply the structure container as an argument so that we can get the structures that we need.
// The connectivity container then needs to a part of the structure container such that we can also
// query the boundary conditions and connectivity.
{
    // useful typedefs
    typedef typename plane_type;
    typedef typename plane_type::length_type    length_type;
    typedef typename position_type;
    typedef typename structure_id_type;
    typedef typename position_structid_pair_type;

    // Note that we assume that the new position is in the plane (dot(pos, unit_z)==0)
    const plane_type origin_plane (sc.get_structure(pos_structure_id.second)->shape());
    boost::array<length_type, 2> half_extends(origin_plane.half_extent());
    const length_type origin_length_x (2*half_extends[0]);
    const length_type origin_length_y (2*half_extends[1]);

    // TODO do cyclic transpose here first
    const position_type pos_vector( subtract(pos_structure_id.first, origin_plane.position()) );
    const length_type component_x ( dot_product(pos_vector, origin_plane.unit_x()) );
    const length_type component_y ( dot_product(pos_vector, origin_plane.unit_y()) );

    // declare the variables that will be written
    structure_id_type new_id;
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
        const position_type neighbor_plane_unit_inline (sc.get_neighbor_info(pos_structure_id.second, 0).second);

        new_id = neighbor_mapping_[pos_structure_id.second][0].first;
        neighbor_plane_par = component_x * origin_plane.unit_x();
        neighbor_plane_inl = (component_y - origin_length_y)* neighbor_plane_unit_inline;
    }
    else if ( (-half_extends[0] < component_x) && ( component_x < half_extends[0]) &&
              (component_y < -half_extends[1] ) )
    {
        // we are at the 'bottom' of the plane (side nr. 1)
        // get the unit vector pointing from the edge between the two planes to the center of the plane
        const position_type neighbor_plane_unit_inline (sc.get_neighbor_info(pos_structure_id.second, 1).second);

        new_id = neighbor_mapping_[pos_structure_id.second][1].first;
        neighbor_plane_par = component_x * origin_plane.unit_x();
        neighbor_plane_inl = (-component_y - origin_length_y)* neighbor_plane_unit_inline;
    }
    else if ( (-half_extends[1] < component_y) && ( component_y < half_extends[1]) &&
              (component_x < -half_extends[0]) )
    {
        // we are at the 'left' of the plane (side nr. 2)
        // get the unit vector pointing from the edge between the two planes to the center of the plane
        const position_type neighbor_plane_unit_inline (sc.get_neighbor_info(pos_structure_id.second, 2).second);

        new_id = neighbor_mapping_[pos_structure_id.second][2].first;
        neighbor_plane_par = component_y * origin_plane.unit_y();
        neighbor_plane_inl = (-component_x - origin_length_x)* neighbor_plane_unit_inline;
    }
    else if ( (-half_extends[1] < component_y) && ( component_y < half_extends[1]) &&
              (half_extends[0] < component_x) )
    {
        // we are at the 'right' of the plane (side nr. 3)
        // this is a unit vector pointing from the edge between the two planes to the center of the plane
        const position_type neighbor_plane_unit_inline (sc.get_neighbor_info(pos_structure_id.second, 3).second);

        new_id = neighbor_mapping_[pos_structure_id.second][3].first;
        neighbor_plane_par = component_y * origin_plane.unit_y();
        neighbor_plane_inl = (component_x - origin_length_x)* neighbor_plane_unit_inline;
    }
    else
    {
        // we are somewhere in the corners and don't know what to do
        // throw exception?
    }

    const position_type new_pos ( add(sc.get_structure(new_id)->position(),
                                      add(neighbor_plane_par, neighbor_plane_inl)));
    // TODO do normal apply_boundary too?
    return std::make_pair(new_id, new_pos);
}

template<>
inline
cyclic_transpose

*/

///// Functions for Cylinders

//// Functions for Boxes


#endif /* STRUCTURE_CONTAINER_HPP */

