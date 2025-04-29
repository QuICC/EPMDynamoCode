/** \file DenseOperator.cpp
 *  \brief Source of the implementation of a general dense operator
 */

// System includes
//

// External includes
//

// Class include
//
#include "Operators/DenseOperator.hpp"

// Project includes
//

namespace EPMDynamo {

   DenseOperator::DenseOperator(const int rows, const int cols)
      : mIsFactorised(false), mOperator(rows, cols)
   {
   }

   DenseOperator::DenseOperator(const int size, const bool isFactorised)
      : mIsFactorised(isFactorised), mOperator(size, size)
   {
   }

   void DenseOperator::solveEquation(Array& vector)
   {
      // This function should never be called
      assert(false);
   }

   void DenseOperator::inverse()
   {
      // This function should never be called
      assert(false);
   }

}
