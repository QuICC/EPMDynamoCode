/** \file EPMTypedefs.hpp
 *  \brief Definition of some useful typedefs used in the whole project
 */

#define EIGEN_YES_I_KNOW_SPARSE_MODULE_IS_NOT_STABLE_YET
#define EIGEN_INITIALIZE_MATRICES_BY_ZERO

#ifndef EPMTYPEDEFS_HPP
#define EPMTYPEDEFS_HPP

// Config includes
//
#include "Config/SmartPointer.h"

// System includes
//
#include <algorithm>
#include <complex>
#include <assert.h>
#include <iostream>

// External includes
//
#include <Eigen/Core>
#include <Eigen/Sparse>

// Project includes
//

namespace EPMDynamo {

   /**
    * @name Basic scalar types typedefs
    */
   //@{
   /// Typedef for integer type value
   typedef int EPMInt;
   /// Typedef for floating point type value
   typedef double EPMFloat;
   /// Typedef for complex type value
   typedef std::complex<EPMFloat> EPMComplex;
   //@}

   /**
    * @name Array/Matrix types typedefs
    */
   //@{
   /// Typedef for a matrix of integer values
   typedef Eigen::Matrix<EPMInt, Eigen::Dynamic, Eigen::Dynamic>   MatrixI;
   /// Typedef for a matrix of float values
   typedef Eigen::Matrix<EPMFloat, Eigen::Dynamic, Eigen::Dynamic>   Matrix;
   /// Typedef for a matrix of complex values
   typedef Eigen::Matrix<EPMComplex, Eigen::Dynamic, Eigen::Dynamic>   MatrixZ;
   /// Typedef for an array of boolean values
   typedef Eigen::Matrix<bool, Eigen::Dynamic, 1>   ArrayB;
   /// Typedef for an array of integer values
   typedef Eigen::Matrix<EPMInt, Eigen::Dynamic, 1>   ArrayI;
   /// Typedef for an array of float values
   typedef Eigen::Matrix<EPMFloat, Eigen::Dynamic, 1>   Array;
   /// Typedef for an array of complex values
   typedef Eigen::Matrix<EPMComplex, Eigen::Dynamic, 1>   ArrayZ;
   //@}
   
   /**
    * @name SparseMatrix typedefs
    */
   //@{
   /// Typedef for a sparse matrix of float values
   typedef Eigen::SparseMatrix<EPMFloat>  SparseMatrix;
   //@}

   /**
    * @name RTP fields/scalars components typedefs
    */
   //@{
   /// Typedef for a matrix of as description for a spherical shell 
   typedef Matrix    SphericalShell;
   /// Typedef for a matrix of as description for a spherical radius 
   typedef Array  SphericalRadius;
   //@}

   /**
    * @name Full Spectral decomposition fields/scalars components typedefs
    */
   //@{
   /// Typedef for a matrix of as description of the set of coefficient for a spectral decomposition at given degree
   typedef MatrixZ  SpectralSHLShell;
   //@}

   /**
    * @name Spherical harmonics decomposition fields/scalars components typedefs
    */
   //@{
   /// Typedef for a array of as description for a spherical shell 
   typedef ArrayZ SHDegrees;
   /// Typedef for a matrix of as description for a spherical shell 
   typedef MatrixZ FDSHLShell;
   //@}
   
   /**
    * @name Array/Matrix Maps typedefs
    */
   //@{
   /// Typedef for a matrix map 
   typedef Eigen::Map<Matrix>   MatrixMap;
   /// Typedef for a complex matrix map 
   typedef Eigen::Map<MatrixZ>   MatrixZMap;
   //@}

   /**
    * @name Boost Smart pointer typedefs
    */
   //@{
   /// Typedef for an smart reference counting pointer on an array of real values
   typedef EPMSHARED_PTR<Array>   SmartArray;
   //@}

   /**
    * @name Timestepping typedefs
    */
   //@{
   /// Typedef for a RHS timestepping matrix
   typedef Matrix  TimestepMatrix;
   //@}
   
   /// Enumb of the different type of predefined timestep controllers
   enum TimestepCtrlTypes {ElementaryCtrl, PI42Ctrl, H211BCtrl};
}

#endif // EPMTYPEDEFS_HPP
