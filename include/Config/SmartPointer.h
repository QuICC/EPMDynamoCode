/** \file SmartPointer.h
 *  \brief Macros/includes depending on the selected smart pointer implementation.
 *
 * Depending on the selected implementation a typedef for the shared_ptr will be provided and the required header files included.
 */

#ifndef SMARTPOINTER_H
#define SMARTPOINTER_H

// Boost Version
#ifdef EPMDYNAMO_BOOST_SMARTPTR
   #include <boost/shared_ptr.hpp>
/**
 * \def EPMSHARED_PTR
 * Macro allowing to use different implementations of the shared_ptr smart pointer.
 * Sets up the Boost version.
 */
   #define EPMSHARED_PTR boost::shared_ptr
#endif //EPMDYNAMO_BOOST_SMARTPTR

// TR1 Version
#ifdef EPMDYNAMO_TR1_SMARTPTR
   #include <tr1/memory>
/**
 * \def EPMSHARED_PTR
 * Macro allowing to use different implementations of the shared_ptr smart pointer.
 * Sets up the TR1 version.
 */
   #define EPMSHARED_PTR std::tr1::shared_ptr
#endif //EPMDYNAMO_TR1_SMARTPTR

// C++0x Version
#ifdef EPMDYNAMO_CXX0X_SMARTPTR
   #include <memory>
/**
 * \def EPMSHARED_PTR
 * Macro allowing to use different implementations of the shared_ptr smart pointer.
 * Sets up the C++0x version.
 */
   #define EPMSHARED_PTR std::shared_ptr
#endif //EPMDYNAMO_CXX0X_SMARTPTR

#endif // SMARTPOINTER_H
