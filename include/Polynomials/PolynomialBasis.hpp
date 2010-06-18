/** \file PolynomialBasis.hpp
 *  \brief Implementation of a general basis of polynomials
 */

#ifndef POLYNOMIALBASIS_HPP
#define POLYNOMIALBASIS_HPP

// System includes
//
#include <vector>

// External includes
//

// Project includes
//
#include "General/EPMTypedefs.hpp"

namespace EPMDynamo {

   /**
    * @brief This class is used for putting together several polynomials to form a basis
    *
    * A polynomialBasis is simply a vector of Polynomials computed a different values of 
    * the fixed parameter. There isn't really anything fancy about this class.
    *
    * \tparam TPoly Type of the polynomial
    *
    * \bug Review this implementation, make it possible to store less polynomials (ie for example odd and even)
    */
   template <typename TPoly> class PolynomialBasis
   {
      public:
         /// Typedef for the basis polynomial type
         typedef TPoly   PolynomialType;

         /**
         * @brief Constructs the polynomial basis from the given values
         *
         * This constructor doesn't do much. It only sets the size and the parameters. 
         * Filling the basis vector will be the responsability of the derived class.
         *
         * @param params   Array of the fixed parameters
         * @param nX       Size of the grid
         * @param nN       Size of the polynomial
         */
         PolynomialBasis(const ArrayI& params, const int nX, const int nN);

         /**
         * @brief Simple empty destructor
         */
         virtual ~PolynomialBasis() {};

         /**
         * @brief Get the size of the basis
         */
         int basisN() const;

         /**
         * @brief Get the number of polynomials per Polynomial object
         */
         int polyN() const;

         /**
         * @brief Size of the grid
         */
         int gridN() const;

         /**
         * @brief Array of parameters of the basis
         *
         * This is the array of values of the fixed parameter of the polynomials
         */
         const ArrayI& parameters() const;

         /**
          * @brief Get polynomial by its position index
          *
          * @param i Position in basis of the requested polynomial
          */
         const TPoly& at(const int i) const ;

      protected:

         /**
          * @brief Holds the set of polynomials defining the basis
          */
         std::vector<TPoly>    mBasis;

      private:
         /**
          * @brief Size of the basis
          */
         int   mBasisN;

         /**
          * @brief Size of the grid
          */
         int   mGridN;

         /**
          * @brief Order of the polynomial
          */
         int   mPolyN;

         /**
          * @brief Holds the set of parameters defining the basis
          */
         ArrayI   mParameters;

   };

   template <typename TPoly> inline int PolynomialBasis<TPoly>::basisN() const
   {
      return this->mBasisN;
   }

   template <typename TPoly> inline int PolynomialBasis<TPoly>::polyN() const
   {
      return this->mPolyN;
   }

   template <typename TPoly> inline int PolynomialBasis<TPoly>::gridN() const
   {
      return this->mGridN;
   }

   template <typename TPoly> inline const ArrayI& PolynomialBasis<TPoly>::parameters() const
   {
      return this->mParameters;
   }

   template <typename TPoly> inline const TPoly& PolynomialBasis<TPoly>::at(const int i) const 
   {
      return this->mBasis.at(i);
   }

   template <typename TPoly> PolynomialBasis<TPoly>::PolynomialBasis(const ArrayI& params, const int nX, const int nN)
      : mBasisN(params.size()), mGridN(nX), mPolyN(nN), mParameters(params)
   {
   }

}

#endif // POLYNOMIALBASIS_HPP
