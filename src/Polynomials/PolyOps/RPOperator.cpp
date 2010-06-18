/** \file RPOperator.cpp
 *  This file contains the definitions of the RPOperator class methods
 */

// System includes
//

// External includes
//

// Class include
//
#include "Polynomials/PolyOps/RPOperator.hpp"

// Project includes
//

namespace EPMDynamo {

   RPOperator::RPOperator(const Array &right, const Matrix &poly)
      : mrPoly(poly), mrR(right)
   {
   }

   RPOperator::~RPOperator()
   {
   }

}
