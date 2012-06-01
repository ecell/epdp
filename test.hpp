// Get_info
// info is:
// -reaction volume
// -product position(s)/structure_id(s)

// In BD, the local situation is checked by determining the overlaps.
// -the reactant particle can overlap with other particles
//  The overlapping particles can be on: -the current structure
//                                       -the parent structure (or the parent of the parent)
//                                       -a peer structure
//                                       -a child structure (or the child of a child)
// -the reactant particle can overlap with other surface.
//  the overlapping surfaces are always: -the current structure (where the particle lives)
//                                       -the parent structure (the parent always encompasses the child surface)
//  the overlapping surface can also be: -a peer structure
//                                       -a child structure

// For the particles we can therefore have three situations. The other particles is:
// -at the same hierarchical level -> need to decide on reaction rule where to place the particle
// -below the current level        -> the structure of the new particle is the structure of 'the other' particle
// -above the current level        -> the structure of the new particle is the current structure (same as previous but particles swapped)


// particles are in the local vicinity
// get the reaction rule for the situation
// check that the reaction rule defines an implemented mechanism/function
// -the product species must be able to live on one of the reactant surfaces ( product_species.STid == particle1.structure.sid/particle2.structure.sid )
// -when calling the right function, the signature for the reaction must exist.
//  (for example: P(region) + P(cylinder) -> P(region) is not supported, whereas P(region) + P(cylinder) -> P(cylinder) is.
// call the right function

if (particles on connected structures of same structure_type)
    // transition/edge_crossing
else
    // particles live on same/different structure_types that are not connected
    determine ordering of structures on which particles live // assume this to be unimportant or unambiguous
    get_product_info -> (position, structure_id)


// Here we assume that the peers/parent conditions are checked in the functions themselves. The relationship between the structure is assumed to be unambiguous or unimportant.
get_pair_product (cuboidal_region_type const& reactant1_region,       cuboidal_region_type const& reactant2_region, cuboidal_region_type const& product_region);            // OK, structures must be peers
get_pair_product (planar_surface_type const& reactant1_plane,         planar_surface_type const& reactant2_plane,   planar_surface_type const& product_plane);              // OK, structures must be peers
get_pair_product (cylindrical_surface_type const& reactant1_cylinder, cylindrical_surface_type const& reactant2_cylinder, cylindrical_surface_type const& product_cylinder);// OK, structures must be peers
get_pair_product (disk_surface_type const& reactant1_disk,            disk_surface_type const& reactant2_disk, ...);                                                        // Illegal


get_pair_product (planar_surface_type const& reactant1_plane,         cuboidal_region_type const& reactant2_region, planar_surface_type const& product_plane);         // Ok, plane must be below region
get_pair_product (cuboidal_region_type const& reactant1_region,       planar_surface_type const& reactant2_plane,   planar_surface_type const& product_plane);         // OK, plane must be below region


get_pair_product (cylindrical_surface_type const& reactant1_cylinder, cuboidal_region_type const& reactant2_region,       cylindrical_surface_type const& product_cylinder);// Ok,cylinder must be below region
get_pair_product (cuboidal_region_type const& reactant1_region,       cylindrical_surface_type const& reactant2_cylinder, cylindrical_surface_type const& product_cylinder);// OK,cylinder must be below region


get_pair_product (disk_surface_type const& reactant1_disk,            cuboidal_region_type const& reactant2_region, disk_surface_type const& product_disk);     // Ok, disk must be below region (Cap/Sink + Bulk -> Cap/Sink)
get_pair_product (cuboidal_region_type const& reactant1_region,       disk_surface_type const& reactant2_disk,      disk_surface_type const& product_disk);     // Ok, disk must be below region


get_pair_product (cylindrical_surface_type const& reactant1_cylinder, planar_surface_type const& reactant2_plane,         ...);         // Illegal, maybe later ok if the structures are peers
get_pair_product (planar_surface_type const& reactant1_plane,         cylindrical_surface_type const& reactant2_cylinder, ...);         // Illegal


get_pair_product (disk_surface_type const& reactant1_disk,            planar_surface_type const& reactant2_plane, ...);                 // Illegal, maybe later ok if the structures are peers
get_pair_product (planar_surface_type const& reactant1_plane,         disk_surface_type const& reactant2_disk,    ...);                 // Illegal


get_pair_product (disk_surface_type const& reactant1_disk,            cylindrical_surface_type const& reactant2_cylinder, disk_surface_type const& product_disk);   // OK, disk must be below cylinder or peer of cylinder
get_pair_product (cylindrical_surface_type const& reactant1_cylinder, disk_surface_type const& reactant2_disk, disk_surface_type const& product_disk);  // OK, disk must be below cylinder or peer of cylinder







// check if the local situation allows the reaction rule to be executed
reactant_structure = get_structure(reactant.structure_id())
if ( product_species.sid() == reactant_structure.parent_structure.sid() )
{
    // The reactant_species lives on the reactant structure
    // Then the reaction rule A(struct_a) -> B(struct_b) where struct_b is the parent of struct_a can be executed (dissociation)
    //  because the parent structure is of the type denoted in the reaction rule.
    get_info_to_parent(reactant_structure, parent_structure)
}
else if ( )
{
    // The reaction rule A(struct_a) -> B(struct_b) where struct_a and struct_b are peers can be executed (transfer).
}
else if ( )
{
    // The reaction A(struct_a) + struct_b -> B(struct_b), where struct_b 
}

// Single_reactions from a structure to another structure
// When the target structure is the parent structure of the source structure
// The structures are the structures that the particle lives on

// These are functions for a single particle that moves from one structure to another
// Note the two structures need to be 'close'. This means that:
// -the target structure is the parent of the source structure (the parent always surrounds the child and is therefore always close)
// -the target structure is the child of the source structure and is 'close'
// -the target structure and source structures are peers and 'close'
// -the two structure are actually the very same structure.

// What to do:
// The two structures are of the same type and are actually the same structure -> no changes (simple decay)
// The two structures are of the same type (but are not the same) -> direct 'transfer' from one structure to the other (hopping or whatever)
// The two structures are of different types -> association/dissociation reaction
// We assume that the relationship between the structures is clear/unambigous or can be determined by the function itself.
get_single_product (cuboidal_region_type const& region,       cuboidal_region_type const& parent_region);       // OK, hopping (Illegal for now) or simple decay (structure are the same one)
get_single_product (planar_surface_type const& plane,         planar_surface_type const& parent_plane);         // OK, hopping (Illegal for now) or simple decay (structure are the same one)
get_single_product (cylindrical_surface_type const& cylinder, cylindrical_surface_type const& parent_cylinder); // OK, hopping (Illegal for now) or simple decay (structure are the same one)
get_single_product (disk_surface_type const& disk,            disk_surface_type const& parent_disk);            // OK, hopping (Illegal for now) or simple decay (structure are the same one)

get_single_product (planar_surface_type const& plane,         cuboidal_region_type const& parent_region);       // OK, dissociation (region must be parent of plane)
get_single_product (cuboidal_region_type const& region,       planar_surface_type const& parent_plane);         // OK, association

get_single_product (cylindrical_surface_type const& cylinder, cuboidal_region_type const& parent_region);       // OK, dissociation (region must be parent of cylinder)
get_single_product (cuboidal_region_type const& region,       cylindrical_surface_type const& parent_cylinder); // OK, association

get_single_product (disk_surface_type const& disk,            cuboidal_region_type const& parent_region);       // OK, dissociation (Cap -> Bulk)
get_single_product (cuboidal_region_type const& region,       disk_surface_type const& parent_disk);            // Illegal for now (association from bulk)

get_single_product (cylindrical_surface_type const& cylinder, planar_surface_type const& parent_plane);         // Illegal
get_single_product (planar_surface_type const& plane,         cylindrical_surface_type const& parent_cylinder); // Illegal
get_single_product (disk_surface_type const& disk,            planar_surface_type const& parent_plane);         // Illegal
get_single_product (planar_surface_type const& plane,         disk_surface_type const& parent_disk);            // Illegal

get_single_product (disk_surface_type const& disk,            cylindrical_surface_type const& parent_cylinder); // OK Sink/Cap -> Cylinder
get_single_product (cylindrical_surface_type const& cylinder, disk_surface_type const& parent_disk);            // OK Cylinder -> Sink/Cap






get_product_pair (cuboidal_region_type const& reactant_region,       cuboidal_region_type const& product1_region, cuboidal_region_type const& product2_region);            // OK, structures must be the same
get_product_pair (planar_surface_type const& reactant_plane,         planar_surface_type const& product1_plane,   planar_surface_type const& product2_plane);              // OK, structures must be the same
get_product_pair (cylindrical_surface_type const& reactant_cylinder, cylindrical_surface_type const& product1_cylinder, cylindrical_surface_type const& product2_cylinder);// OK, structures must be the same
get_product_pair (disk_surface_type const& reactant_disk,            disk_surface_type const& product1_disk, disk_surface_type const& product2_disk);                      // Illegal
                
                
get_product_pair (planar_surface_type const& reactant_plane,         cuboidal_region_type const& product1_region, planar_surface_type const& product_plane);         // Ok, plane must be below region
get_product_pair (cuboidal_region_type const& reactant_region,       planar_surface_type const& product1_plane,   planar_surface_type const& product_plane);         // OK, plane must be below region
                
                
get_product_pair (cylindrical_surface_type const& reactant_cylinder, cuboidal_region_type const& reactant2_region,       cylindrical_surface_type const& product_cylinder);// Ok,cylinder must be below region
get_product_pair (cuboidal_region_type const& reactant_region,       cylindrical_surface_type const& reactant2_cylinder, cylindrical_surface_type const& product_cylinder);// OK,cylinder must be below region
                
                
get_product_pair (disk_surface_type const& reactant_disk,            cuboidal_region_type const& reactant2_region, disk_surface_type const& product_disk);     // Ok, disk must be below region (Cap/Sink + Bulk -> Cap/Sink)
get_product_pair (cuboidal_region_type const& reactant_region,       disk_surface_type const& reactant2_disk,      disk_surface_type const& product_disk);     // Ok, disk must be below region
                
                
get_product_pair (cylindrical_surface_type const& reactant_cylinder, planar_surface_type const& reactant2_plane,         ...);         // Illegal, maybe later ok if the structures are peers
get_product_pair (planar_surface_type const& reactant_plane,         cylindrical_surface_type const& reactant2_cylinder, ...);         // Illegal
                
                
get_product_pair (disk_surface_type const& reactant_disk,            planar_surface_type const& reactant2_plane, ...);                 // Illegal, maybe later ok if the structures are peers
get_product_pair (planar_surface_type const& reactant_plane,         disk_surface_type const& reactant2_disk,    ...);                 // Illegal
                
                
get_product_pair (disk_surface_type const& reactant_disk,            cylindrical_surface_type const& reactant2_cylinder, disk_surface_type const& product_disk);   // OK, disk must be below cylinder or peer of cylinder
get_product_pair (cylindrical_surface_type const& reactant_cylinder, disk_surface_type const& reactant2_disk, disk_surface_type const& product_disk);  // OK, disk must be below cylinder or peer of cylinder



class Structure
{
    virtual position_type method(structure_type const& target_structure, position_type const& position) = 0;
}

class PlanarSurface : Structure
{

    virtual position_type method(structure_type const& target_structure, position_type const& position)
    {
        return target_structure.get_pos(*this, position);
    }

    template <Tstruct_>
    virtual position_type get_pos(Tstruct_ const& origin_structure, position_type const& position)
    {
        return ::get_pos(origin_structure, *this, position);
    }

}

class CylindricalSurface
{

    virtual position_type method(structure_type const& target_structure, position_type const& position)
    {
        return target_structure.get_pos(*this, position);
    }

    template <Tstruct_>
    virtual position_type get_pos(Tstruct_ const& origin_structure, position_type const& position)
    {
        return ::get_pos(origin_structure, *this, position);
    }

}

#include "PlanarSurface.hpp"
#include "CylindricalSurface.hpp"
#include "DiskSurface.hpp"

template <Ttraits_>
position_type get_pos(PlanarSurface<Ttraits_> const& origin_structure,
                      PlanarSurface<Ttraits_> const& target_structure, typename Ttraits_::position_type const& old_pos)
{
    typedef typename Ttraits_::structure_id_type structure_id_type;
    typedef typename Ttraits_::position_type     position_type;

    const structure_id_type id();
    return std::make_pair(target_structure.project_point(old_pos), target_structure.id());
}

get_pos(cylinder_surface_type, cylinder_surface_type)
{
    throw exception;
}
template <Ttraits_>
inline PlanarSurface<Ttraits_>::position_structid_pair_type
get_pos(CylindricalSurface<Ttraits_> origin_structure,
        planar_surface_type);


get_pos(planar_surface_type, cylinder_surface_type);


newBDthing
    new_pos = source_structure.method(target_structure, old_pos);

