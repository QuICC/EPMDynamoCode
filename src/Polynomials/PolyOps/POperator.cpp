/** \file POperator.cpp
 *  This file contains the definitions of the POperator class methods
 */

// System includes
//

// External includes
//

// Class include
//
#include "Polynomials/PolyOps/POperator.hpp"

// Project includes
//

namespace EPMDynamo {

   POperator::POperator(const Matrix &poly)
      : mrPoly(poly)
   {
   }

   POperator::~POperator()
   {
   }

}
