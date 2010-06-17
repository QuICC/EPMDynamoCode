/** \file RadialTransform.hpp
 *  \brief Implementation of a general radial transform
 */

#ifndef RADIALTRANSFORM_HPP
#define RADIALTRANSFORM_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "Domain/Truncation.hpp"
#include "Polynomials/RadialBasis.hpp"
#include "Transforms/Polynomial/QuadratureMultiplier.hpp"

namespace EPMDynamo {

   /**
    * \brief Implementation of a general radial transform
    *
    * \tparam TPolynomial Radial polynomial basis
    */
   template <typename TPolynomial> class RadialTransform: public QuadratureMultiplier<RadialBasis<TPolynomial> >
   {
      public:
         /**
          * @brief Constructs the transform object
          *
          * This object simply holds a polynomials basis
          * and the methods required to perform forward and backward transforms
          *
          * @param pTrunc Truncation information
          */
         RadialTransform(SmartTruncation pTrunc);

         /**
          * @brief Simple empty destructor
          */
         virtual ~RadialTransform() {};

      protected:

      private:
   };

   template <typename TPolynomial> RadialTransform<TPolynomial>::RadialTransform(SmartTruncation pTrunc)
      : QuadratureMultiplier<RadialBasis<TPolynomial> >(pTrunc->local()->spec()->lArray(), pTrunc->sim()->rad()->nR(), pTrunc->sim()->rad()->nN())
   {
   }
}

#endif // RADIALTRANSFORM_HPP
