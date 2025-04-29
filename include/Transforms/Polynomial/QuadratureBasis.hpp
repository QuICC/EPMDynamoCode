/** \file QuadratureBasis.hpp
 *  \brief Base of the implementation of a quadrature transform 
 */

#ifndef QUADRATUREBASIS_HPP
#define QUADRATUREBASIS_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "General/EPMTypedefs.hpp"

namespace EPMDynamo {

   /**
    * \brief Base of the implementation of a quadrature transform 
    *
    * \tparam TBasis Basis used in quadrature
    */
   template<typename TBasis>  class QuadratureBasis
   {
      public:
         /**
         * @brief Constructs the basis required for the quadrature
         *
         * This templated class is only reponsible for initialisating
         * the basis
         *
         * @param params   Array of the fixed parameters
         * @param nX       Size of the grid
         * @param nN       Size of the polynomial
         */
         QuadratureBasis(const ArrayI& params, const int nX, const int nN);

         /**
         * @brief Simple empty destructor
         */
         virtual ~QuadratureBasis() {};

         /**
          * @brief Set the basis
          */
         TBasis&  basis();

         /**
          * @brief Get the basis
          */
         const TBasis&  basis() const;

      protected:
         /**
          * @brief Size of the basis
          */
         int mN;

         /**
          * @brief Pointer to quadrature basis
          */
         TBasis   mBasis;
      
      private:
   };

   template<typename TBasis>  QuadratureBasis<TBasis>::QuadratureBasis(const ArrayI& params, const int nX, const int nN)
      : mN(params.size()), mBasis(params, nX, nN)
   {
   }

   template<typename TBasis> TBasis& QuadratureBasis<TBasis>::basis()
   {
      return this->mBasis;
   }

   template<typename TBasis> const TBasis& QuadratureBasis<TBasis>::basis() const
   {
      return this->mBasis;
   }

}

#endif // QUADRATUREBASIS_HPP
