#ifndef BINDING_MODEL_HPP 
#define BINDING_MODEL_HPP

#include <boost/python.hpp>
#include "peer/utils.hpp"

namespace binding {



template<typename Tmodel_>
struct rfr_structure
{
   typedef Tmodel_ impl_type;
   typedef typename boost::remove_reference<typename impl_type::structure_types_range>::type range_type;
   typedef typename boost::range_const_iterator<range_type>::type result_type;

   typedef peer::wrappers::stl_iterator_wrapper<result_type, boost::python::object> wrapper_type;

   static PyObject* create_wrapper(boost::python::back_reference<const impl_type &> backref)
   {
      wrapper_type::__class_init__(typeid(result_type).name());
      return wrapper_type::create(backref.get().get_structure_types(), backref.source());
   }
};



////// Registering master function
template<typename Tparticle_model_>
inline void register_particle_model_class(char const* name)
{
    using namespace boost::python;
    typedef Tparticle_model_ impl_type;

    // defining the python class
    class_<impl_type, bases<typename impl_type::base_type>, boost::noncopyable>(name)
        .def("add_structure_type", &impl_type::add_structure_type)
        .def("get_structure_type_by_id", &impl_type::get_structure_type_by_id)
        .def("get_def_structure_type_id", &impl_type::get_def_structure_type_id)
        .add_property("structure_types",
           boost::python::make_function(&rfr_structure<impl_type>::create_wrapper))
       ;
}

} // namespace binding

#endif /* BINDING_MODEL_HPP */
