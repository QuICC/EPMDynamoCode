/** \file TorPolRadialOperatorBase.hpp
 *  \brief Base of the implementation of a Toroidal/Poloidal set of radial operators
 */

#ifndef TORPOLRADIALOPERATORBASE_HPP
#define TORPOLRADIALOPERATORBASE_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "General/EPMTypedefs.hpp"
#include "Polynomials/RadialOperatorBase.hpp"

namespace EPMDynamo {

   /**
    * \brief Base of the implementation of a Toroidal/Poloidal set of radial operators
    *
    * \tparam TPolynomial Type of the underlying polynomial
    */
   template <typename TPolynomial> class TorPolRadialOperatorBase: public RadialOperatorBase<TPolynomial>
   {
      public:
         /**
          * @brief Constructs the polynomials and sets up addition operators
          *
          * This constructor computes the values of the grid and the weights and then
          * compute the polynomial values based on these
          *
          * @param pL      Polynomial L order
          * @param nX      Size of the grid
          * @param nN      Size of polynomial basis
          */
         TorPolRadialOperatorBase(const int pL, const int nX, const int nN);

         /**
          * @brief Constructs the polynomials with given points and weighs and sets up
          * addition operators
          *
          * The values of the polynomial are computed on the given grid and weighted with
          * the given weights.
          *
          * @param pL      Polynomial L order
          * @param grid    Provided grid array
          * @param nN      Size of polynomial basis
          * @param weights Provided weights array
          */
         TorPolRadialOperatorBase(const int pL, const SmartArray grid, const int nN, const SmartArray weights);

         /**
          * @brief Simple empty destructor
          */
         virtual ~TorPolRadialOperatorBase() {};

         /**
          * @brief Get the energy matrix for a poloidal component
          */
         const Matrix&  polEWeights() const;

      protected:
         /**
          * @brief Store the QST S component to curl curl projection integrator 
          */
         Matrix   mS2CurlCurlIntegrator;

         /**
          * @brief Store the Energy integral weights of the polynomials for the poloidal
          *    component
          */
         Matrix    mPolEWeights;

      private:
         /**
          * @brief Compute the integrators for the toroidal/poloidal
          *          projection of the non linear terms
          */
         void computeNTermsIntegrators();

         /**
          * @brief Compute the energy integral factors
          *
          * This computes and stores the values \f$ \int_{0}^{1} r^2 P d r\f$
          */
         void computeEnergyIntegrals();
   };

   //
   // Inline definitions of the spectral operators
   //

   template <typename TPolynomial> inline const Matrix& TorPolRadialOperatorBase<TPolynomial>::polEWeights() const
   {
      return this->mPolEWeights;
   }

   template <typename TPolynomial> TorPolRadialOperatorBase<TPolynomial>::TorPolRadialOperatorBase(const int pL, const int nX, const int nN)
      : RadialOperatorBase<TPolynomial>(pL, nX, nN), mS2CurlCurlIntegrator(nN, nX), mPolEWeights(nN, nN)
   {
      // Compute the NTerms integrators
      this->computeNTermsIntegrators();

      // Compute the energy integrals
      this->computeEnergyIntegrals();
   }

   template <typename TPolynomial> TorPolRadialOperatorBase<TPolynomial>::TorPolRadialOperatorBase(const int pL, const SmartArray grid, const int nN, const SmartArray weights)
      : RadialOperatorBase<TPolynomial>(pL, grid, nN, weights), mS2CurlCurlIntegrator(nN, this->gridN()), mPolEWeights(nN, nN)
   {
      // Compute the NTerms integrators
      this->computeNTermsIntegrators();

      // Compute the energy integrals
      this->computeEnergyIntegrals();
   }

   template <typename TPolynomial> void TorPolRadialOperatorBase<TPolynomial>::computeNTermsIntegrators()
   {
      this->mS2CurlCurlIntegrator = -this->sll_1()*(this->grid().asDiagonal()*this->wPoly()*this->diff(1)*this->mR_1.asDiagonal()*this->wPoly());
   }

   template <typename TPolynomial> void TorPolRadialOperatorBase<TPolynomial>::computeEnergyIntegrals()
   {
      int ptsN = 2*this->polyN() + this->l() + 1;

      // Initialise the grid generator
      LegendreGrid  generator(ptsN);

      // Storage for the temporary grid
      SmartArray eGrid(new Array(ptsN));
      SmartArray eWeights(new Array(ptsN));

      // get computed grid
      (*eGrid) = generator.grid();

      // get computed weights
      (*eWeights) = generator.weights();

      // Adapt weights to stretched interval
      *eWeights = (*eWeights)/2.0;

      // Strech grid points to 0, 1 interval
      eGrid->array() += 1.0;
      *eGrid = (*eGrid)/2.0;

      TPolynomial    tmpPoly(this->l(), eGrid, this->polyN(), eWeights); 

      Matrix   tmp(this->polyN(), ptsN);

      EPMFloat lfactor = static_cast<EPMFloat>(this->l()*(this->l()+1));
      EPMFloat correctNorm;

      for(int n=0; n < this->polyN(); ++n)
      {
         for(int i=0; i < this->polyN(); ++i)
         {
            correctNorm = this->normalisation()(i)*this->normalisation()(n)/(tmpPoly.normalisation()(i)*tmpPoly.normalisation()(n));

            this->mPolEWeights(i,n) = (lfactor*(lfactor+1.0))*(tmpPoly.poly().row(n).array()*tmpPoly.poly().row(i).array()).matrix().dot (*eWeights);
            this->mPolEWeights(i,n) += lfactor*2.0*(tmpPoly.poly().row(n).array()*tmpPoly.diff(1).row(i).array()).matrix().dot ((eGrid->array()*eWeights->array()).matrix());
            this->mPolEWeights(i,n) += lfactor*(tmpPoly.diff(1).row(n).array()*tmpPoly.diff(1).row(i).array()).matrix().dot ((eGrid->array().pow(2)*eWeights->array()).matrix());

            this->mPolEWeights(i,n) *= correctNorm;
         }  
      }
   }

}

#endif // TORPOLRADIALOPERATORBASE_HPP
