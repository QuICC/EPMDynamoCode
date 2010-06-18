/** \file QuadratureRule.hpp
 *  \brief Implementation of the computation of grid points and weights for a quadrature formul
 */

#ifndef QUADRATURERULE_HPP
#define QUADRATURERULE_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "General/EPMTypedefs.hpp"

namespace EPMDynamo {

   /**
    * @brief This class provides the methods to compute a quadrature rule
    *
    * This class provides all the required method to compute quadrature points and 
    * quadrature weights. The linear algebra computations are currently done with lapack
    * but might be done directly with Eigen at a later point. Speed isn't important here 
    * as it is only used in an initialisation step. The only important thing is maximal
    * precision.
    *
    * Documentation on how these computations are done can be found in Golub's 1969 paper
    */
   class QuadratureRule
   {
      public:
         /**
          * @brief Construct the quadrature rule tridiagonal matrix
          *
          * This constructor currently requires the diagonal and subdiagonal arrays of 
          * the matrix as input. It might be rewritten in a simple fashion later on. For
          * exmaple to only require pointers to the recurrence coefficient methods or
          * something similar.
          *
          * @param diagonal      The diagonal of the symmetric tridiagonal eigensystem
          * @param subdiagonal   The subdiagonal of the symmetric tridiagonal eigensystem
          * @param intWeight     The integral weight of the weight function
          */
         QuadratureRule(const Array& diagonal, const Array& subdiagonal, DynamoFloat intWeight = 0.0);

         /**
          * @brief Simple empty destructor
          */
         virtual ~QuadratureRule() {};

         /**
          * @brief Extract grid point from the matrix
          *
          * @param grid    Array to store grid points
          */
         void computeGrid(Array& grid);

         /**
          * @brief Extract grid points and weights from the matrix
          *
          * @param grid       Array to store grid points
          * @param weights    Array to store weights
          */
         void computeGridAndWeights(Array& grid, Array& weights);
         
      protected:

      private:
         /**
          * @brief Diagonal of the eigensystem
          */
         Array    mDiagonal;

         /**
          * @brief Subdiagonal of the eigensystem
          */
         Array    mSubdiagonal;

         /**
          * @brief Integral weight of the weight function
          */
         DynamoFloat    mIntWeight;
   };

}

#endif // QUADRATURERULE_HPP
