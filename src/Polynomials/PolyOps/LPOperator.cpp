/** \file LPOperator.cpp
 *  This file contains the definitions of the LPOperator class methods
 */

// System includes
//

// External includes
//

// Class include
//
#include "Polynomials/PolyOps/LPOperator.hpp"

// Project includes
//

namespace EPMDynamo {

   LPOperator::LPOperator(const Array &left, const Matrix &poly)
      : mrPoly(poly), mrL(left)
   {
   }

   LPOperator::~LPOperator()
   {
   }

}
