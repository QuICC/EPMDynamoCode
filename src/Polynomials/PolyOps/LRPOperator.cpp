/** \file LRPOperator.cpp
 *  This file contains the definitions of the LRPOperator class methods
 */

// System includes
//

// External includes
//

// Class include
//
#include "Polynomials/PolyOps/LRPOperator.hpp"

// Project includes
//

namespace EPMDynamo {

   LRPOperator::LRPOperator(const Array &left, const Array &right, const Matrix &poly)
      : mrL(left), mrR(right), mrPoly(poly)
   {
   }

   LRPOperator::~LRPOperator()
   {
   }

}
