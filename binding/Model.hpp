#ifndef BINDING_MODEL_HPP 
#define BINDING_MODEL_HPP

#include <boost/python.hpp>
#include <boost/python/docstring_options.hpp>
#include "peer/utils.hpp"

namespace binding {

template<typename Timpl>
static void Model___setitem__(Timpl* model, std::string const& key, std::string const& value)
{
    (*model)[key] = value;
}


template<typename Tmodel_>
struct rfr_attr
{
   typedef Tmodel_ impl_type;
   typedef typename boost::remove_reference<typename impl_type::attributes_range>::type range_type;
   typedef typename boost::range_const_iterator<range_type>::type result_type;

   typedef peer::wrappers::stl_iterator_wrapper<result_type, boost::python::object> wrapper_type;

   static PyObject* create_wrapper(boost::python::back_reference<const impl_type &> backref)
   {
      wrapper_type::__class_init__(typeid(result_type).name());
      return wrapper_type::create(backref.get().attributes(), backref.source());
   }
};

template<typename Tmodel_>
struct rfr_species
{
   typedef Tmodel_ impl_type;
   typedef typename boost::remove_reference<typename impl_type::species_type_range>::type range_type;
   typedef typename boost::range_const_iterator<range_type>::type result_type;

   typedef peer::wrappers::stl_iterator_wrapper<result_type, boost::python::object> wrapper_type;

   static PyObject* create_wrapper(boost::python::back_reference<const impl_type &> backref)
   {
      wrapper_type::__class_init__(typeid(result_type).name());
      return wrapper_type::create(backref.get().get_species_types(), backref.source());
   }
};



////// Registering master function
template<typename Tmodel_>
inline void register_model_class(char const* name)
{
    using namespace boost::python;
    typedef Tmodel_ impl_type;

    docstring_options doc_options;
    doc_options.disable_signatures();

    // defining the python class
    class_<impl_type, boost::noncopyable>(name)
        .add_property("network_rules",
            make_function(&impl_type::network_rules,
                return_value_policy<reference_existing_object>()))
        .def("add_species_type", &impl_type::add_species_type)
        .def("get_species_type_by_id", &impl_type::get_species_type_by_id)
        .def("__getitem__", (std::string const&(impl_type::*)(std::string const&) const)
                &impl_type::operator[], return_value_policy<copy_const_reference>())
        .def("__setitem__", &Model___setitem__<Tmodel_>)
        .add_property("attributes",
               boost::python::make_function(&rfr_attr<impl_type>::create_wrapper))
        .add_property("species_types",
               boost::python::make_function(&rfr_species<impl_type>::create_wrapper))
       ;
}

} // namespace binding

#endif /* BINDING_MODEL_HPP */

