/** \file AssociatedLegendreOperatorBase.hpp
 *  \brief Base of the implementation of a set of operators based on the associated Legendre polynomial
 */

#ifndef ASSOCIATEDLEGENDREOPERATORBASE_HPP
#define ASSOCIATEDLEGENDREOPERATORBASE_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "Polynomials/AssociatedLegendrePolynomial.hpp"
#include "General/MathConstants.hpp"

namespace EPMDynamo {

   /**
    * \brief Base of the implementation of a set of operators based on the associated Legendre polynomial
    */
   class AssociatedLegendreOperatorBase: public AssociatedLegendrePolynomial
   {
      public:

         /**
          * @brief Constructs the polynomial and the operators from scratch.
          *
          * This constructors computes the values of the grid and the weights and then
          * computes the polynomial values based on these
          *
          * @param nX      Size of the grid
          * @param nL      Number of l modes
          * @param pM      Polynomial m degree
          */
         AssociatedLegendreOperatorBase(const int pM, const int nX, const int nL);

         /**
          * @brief Constructs the polynomial and the operators based on the given grid points and weights
          *
          * The values of the polynomial are computed on the given grid and weighted with
          * the given weights.
          *
          * @param grid    Provided grid array
          * @param nL      Number of l modes
          * @param weights Provided weights array
          * @param pM      Polynomial m degree
          */
         AssociatedLegendreOperatorBase(const int pM, const SmartArray grid, const int nL, const SmartArray weights);

         /**
          * @brief Simple empty destructor
          */
         virtual ~AssociatedLegendreOperatorBase() {};

         /**
          * @brief Get the array of \f$\frac{1}{\sqrt{l(l+1)}}\f$ factors
          */
         const Array& sll_1() const;
         
      protected:
         /**
          * @brief Array holding the 1.0/sin(theta) values
          */
         Array    mSin_1Theta;

         /**
          * @brief Array holding the 1.0/sqrt(l(l+1)) values
          */
         Array    mSll_1;

         /**
          * @brief Array of ones (only needed for CSCS ouput)
          */
         Array    mSin_1Ones;

         /**
          * @brief Special computation for P/Sin if the poles are required (required for CSCS ouput)
          */
         Matrix    mSin_1Poly;

      private:
         /**
          * @brief Compute the 1.0/sin(theta)
          */
         void computeSin_1Theta();

         /**
          * @brief Compute the 1.0/sqrt(l(l+1))
          */
         void computeSll_1();

         /**
          * @brief Compute special operators (for example CSCS stuff)
          *
          * \epmBug Commented out code, needs review and corrections
          */
         void computeSpecialOperators();

   };

   inline const Array&  AssociatedLegendreOperatorBase::sll_1() const
   {
      return this->mSll_1;
   }

}

#endif // ASSOCIATEDLEGENDREOPERATORBASE_HPP
