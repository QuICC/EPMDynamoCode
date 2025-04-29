/** \file RadialBasis.hpp
 *  \brief Implementation of a general radial polynomial basis
 */

#ifndef RADIALBASIS_HPP
#define RADIALBASIS_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "Domain/Truncation.hpp"
#include "Polynomials/PolynomialBasis.hpp"
#include "Polynomials/TorPolRadialOperator.hpp"

namespace EPMDynamo {

   /**
    * @brief This class defines a radial basis
    *
    * This is an implementation of the PolynomialBasis to store a set of 
    * radial operators based on a templated polynomial. There isn't much going
    * on in this class, in mainly serves as a container for the polynomials.
    *
    * \tparam TPolynomial Type of the radial polynomial
    * \tparam TOperatorType Type of the radial operators
    *
    * \epmTodo Review the implementation. As for the PolynomialBasis more flexibilty might be required
    */
   template <typename TPolynomial, template <typename> class TOperatorType = TorPolRadialOperator> class RadialBasis : public PolynomialBasis<TOperatorType<TPolynomial> >
   {
      public:
         /**
         * @brief Constructs the polynomial basis from the given parameters
         *
         * @param parameters  Array of the fixed parameter
         * @param nX          Size of the grid
         * @param nN          Size of the polynomial
         */
         RadialBasis(const ArrayI &parameters, const int nX, const int nN);

         /**
         * @brief Constructs the polynomial based on the given truncation
         *
         * In this case the parameters m are extract from the truncation object.
         *
         * @param pTrunc   Truncation information
         */
         RadialBasis(const SmartTruncation pTrunc);

         /**
         * @brief Simple empty destructor
         */
         virtual ~RadialBasis() {};
         
      protected:

      private:
         /**
          * @brief Fill the basis vector
          */
         void initBasis();
   };

   template <typename TPolynomial, template <typename> class TOperatorType> RadialBasis<TPolynomial, TOperatorType>::RadialBasis(const ArrayI &parameters, const int nX, const int nN)
      : PolynomialBasis<TOperatorType<TPolynomial> >(parameters, nX, nN)
   {
      // Initialise the vector of basis polynomials
      this->initBasis();
   }

   template <typename TPolynomial, template <typename> class TOperatorType> RadialBasis<TPolynomial, TOperatorType>::RadialBasis(const SmartTruncation pTrunc)
      : PolynomialBasis<TOperatorType<TPolynomial> >(pTrunc->local()->spec()->lArray(), pTrunc->sim()->rad()->nR(), pTrunc->sim()->rad()->nN())
   {
      // Initialise the vector of basis polynomials
      this->initBasis();
   }

   template <typename TPolynomial, template <typename> class TOperatorType> void RadialBasis<TPolynomial, TOperatorType>::initBasis()
   {
      // Reserve the space in the vector to avoid multiple reallocations
      this->mBasis.reserve(this->basisN());

      // Full initialisation of the first polynomial (including grid and weights computation)
      this->mBasis.push_back(TOperatorType<TPolynomial>(this->parameters()(0), this->gridN(), this->polyN()));

      //Get pointer to grid points array
      SmartArray grid = this->at(0).pGrid();
      //Get pointer to weights array
      SmartArray weights = this->at(0).pWeights();

      // Initialise the remaining polynomials using the grid and weights pointers
      // In this way the grid and weights are only stored once. The other polynomials
      // simply share the pointer
      for(int i=1; i < this->basisN(); ++i)
      {
         this->mBasis.push_back(TOperatorType<TPolynomial>(this->parameters()(i), grid, this->polyN(), weights));
      }
   }

}

#endif // RADIALBASIS_HPP
