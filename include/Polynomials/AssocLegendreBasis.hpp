/** \file AssocLegendreBasis.hpp
 *  \brief Implementation of a basis of Associated Legendre polynomials
 */

#ifndef ASSOCLEGENDREBASIS_HPP
#define ASSOCLEGENDREBASIS_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "Domain/Truncation.hpp"
#include "Polynomials/PolynomialBasis.hpp"
#include "Polynomials/AssociatedLegendreOperator.hpp"

namespace EPMDynamo {

   /**
    * @brief This class defines a Associated Legendre polynomials basis
    *
    * This is an implementation of the PolynomialBasis to store a set of 
    * Associated Legendre functions. There isn't much going on in this
    * class, in mainly serves as a container for the polynomials.
    */
   class AssocLegendreBasis : public PolynomialBasis<AssociatedLegendreOperator>
   {
      public:
         /**
          * @brief Constructs the polynomial basis from the given parameters
          *
          * @param parameters  Array of the fixed parameter m
          * @param nX          Size of the grid
          * @param nL          Size of the polynomial
          */
         AssocLegendreBasis(const ArrayI &parameters, const int nX, const int nL);

         /**
          * @brief Constructs the polynomial based on the given truncation
          *
          * In this case the parameters m are extract from the truncation object.
          *
          * @param pTrunc      Truncation (will be used to generate parameters)
          */
         AssocLegendreBasis(const SmartTruncation pTrunc);

         /**
          * @brief Simple empty destructor
          */
         virtual ~AssocLegendreBasis() {}:
         
      protected:

      private:

         /**
          * @brief Fill the basis vector
          */
         void initBasis();
   };

}

#endif // ASSOCLEGENDREBASIS_HPP
