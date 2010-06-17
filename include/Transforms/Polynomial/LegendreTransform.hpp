/** \file LegendreTransform.hpp
 *  \brief Implementation of the Legendre transform
 */

#ifndef LEGENDRETRANSFORM_HPP
#define LEGENDRETRANSFORM_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "Transforms/Polynomial/QuadratureMultiplier.hpp"
#include "Domain/Truncation.hpp"
#include "Polynomials/AssocLegendreBasis.hpp"

namespace EPMDynamo {

   /**
    * @brief This class is reponsible for the Legendre Quadrature transforms
    *
    * This is a simple implementation of the QuadratureMultiplier class based
    * on a Associated Legendre polynomial basis.
    */
   class LegendreTransform: public QuadratureMultiplier<AssocLegendreBasis>
   {
      public:
         /**
          * @brief Constructs the transform object
          *
          * This object simply holds a Associated legendre functions basis
          * and the methods required to perform forward and backward transforms
          */
         LegendreTransform(SmartTruncation pTrunc);

         /**
          * @brief Simple empty destructor
          */
         ~LegendreTransform() {};
   };
}

#endif // LEGENDRETRANSFORM_HPP
