

template
// overloaded typedef so as to fix the length of the vectors for different shape types
// -Cylinder n=2
// -Plane n=4
// -Box n=6



template<typename StructureType_>
class ConnectivityContainer
// implements the way structures are connected.
// It allows to to check what structure you enter into when you leave the current one.
{
    typedef typename StructureType_ structure_type;

public:
    position_structid_pair_type apply_boundary (position_structid_pair_type const& pos_structure_id) const
    {
        ::to_internal(structure_type::base_type.shape(), pos);
    }

private:
    // mapping StructureID -> vector<(StructureID, index), max_index+1>

};



class StructureContainer
{

public:
    add_structure (CylindricalSurface)  // can I make this into a template function?
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

    position_structid_pair_type apply_boundary(position_structid_pair_type const& pos_structure_id) const
    // Propagate these to the ParticleContainer such that world.apply_boundary( (pos, struct_id)) can be called
    {
    }
    position_structid_pair_type cyclic_transpose(position_structid_pair_type const& pos_structure_id) const
    // This is essentially the deflect_back function
    {
    }

    remove_structure (StructureID)
    {
        // find StructureID in all ConnectivityContainers -> remove references
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

private:
    structures_range get_structures_range() const
    {
        return structures_range(
        structure_iterator(structure_map_.begin(), structures_second_selector_type()),
        structure_iterator(structure_map_.end(),   structures_second_selector_type()),
        structure_map_.size());
    }
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

protected:
    structure_map                       structure_map_;     // mapping: structure_id -> structure
    per_structure_substructure_id_set   structure_substructures_map_;

    CylindricalSurfaceConnectivityContainer     // This contains a StructureID -> vector <(StructureID, index), 2> map
    PlanarSurfaceConnectivityContainer          // This contains a StructureID -> vector <(StructureID, index), 4> map
    CuboidalReginConnectivityContainer          // This contains a StructureID -> vector <(StructureID, index), 6> map

};
