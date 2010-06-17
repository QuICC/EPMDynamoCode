/** \file LegendreTransform.cpp
 *  \brief Source of implementation of the Legendre transform
 */

// System includes
//

// External includes
//

// Class include
//
#include "Transforms/Polynomial/LegendreTransform.hpp"

// Project includes
//

namespace EPMDynamo {

   LegendreTransform::LegendreTransform(SmartTruncation pTrunc)
      : QuadratureMultiplier<AssocLegendreBasis>(pTrunc->local()->fdsh()->mArray(), pTrunc->sim()->hoz()->nTh(), pTrunc->sim()->hoz()->nL())
   {
   }

}
