#ifndef BINDING_DISK_SURFACE_HPP
#define BINDING_DISK_SURFACE_HPP

#include <boost/python.hpp>

namespace binding {


////// Registering master function
template<typename Timpl>
inline boost::python::objects::class_base register_disk_surface_class(char const *name)
{
    using namespace boost::python;
    typedef Timpl impl_type;

    return class_<impl_type, bases<typename impl_type::base_type::base_type>,
           boost::shared_ptr<impl_type>, boost::noncopyable>(
            name, init<typename impl_type::structure_name_type,
                       typename impl_type::structure_type_id_type,
                       typename impl_type::structure_id_type,
                       typename impl_type::shape_type>())
        .add_property("shape",
            make_function((typename impl_type::shape_type const&(impl_type::*)()const)&impl_type::shape, return_value_policy<return_by_value>()))
        .def("treat_as_barrier", &impl_type::treat_as_barrier)
        .def("treat_as_sink",    &impl_type::treat_as_sink)
        .def("allow_radial_dissociation",  &impl_type::allow_radial_dissociation)
        .def("forbid_radial_dissociation", &impl_type::forbid_radial_dissociation)
        .def("is_barrier",  &impl_type::is_barrier,  return_value_policy<return_by_value>())
        .def("dissociates_radially", &impl_type::dissociates_radially, return_value_policy<return_by_value>())
        ;
}

} // namespace binding

#endif /* BINDING_DISK_SURFACE_HPP */
