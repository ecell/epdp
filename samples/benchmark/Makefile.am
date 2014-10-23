noinst_PROGRAMS = hardbody

AM_CXXFLAGS = -I$(top_srcdir) @BOOST_CPPFLAGS@ @GSL_CFLAGS@ $(PYTHON_INCLUDES)

hardbody_SOURCES = hardbody.cpp ../../NetworkRules.cpp ../../BasicNetworkRulesImpl.cpp ../../Logger.cpp ../../freeFunctions.cpp

hardbody_LDADD = $(GSL_LIBS)
