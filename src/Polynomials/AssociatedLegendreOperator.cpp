/** \file AssociatedLegendreOperator.cpp
 *  \brief Source of the implementation of a set of operators on the associated Legendre poylnomial
 */

// System includes
//

// External includes
//

// Class include
//
#include "Polynomials/AssociatedLegendreOperator.hpp"

// Project includes
//
#include "Simulations/SimulationConstants.hpp"

namespace EPMDynamo {

   AssociatedLegendreOperator::AssociatedLegendreOperator(const int pM, const int nX, const int nL)
      : AssociatedLegendreOperatorBase(pM, nX, nL)
   {
      // Initialise Polynomial operators
      initPOperators();
   }

   AssociatedLegendreOperator::AssociatedLegendreOperator(const int pM, const SmartArray grid, const int nL, const SmartArray weights)
      : AssociatedLegendreOperatorBase(pM, grid, nL, weights)
   {
      // Initialise Polynomial operators
      initPOperators();
   }

   AssociatedLegendreOperator::AssociatedLegendreOperator(const AssociatedLegendreOperator &other)
      : AssociatedLegendreOperatorBase(other)
   {
      // Initialise Polynomial operators
      initPOperators();
   }

   AssociatedLegendreOperator& AssociatedLegendreOperator::operator=(const AssociatedLegendreOperator &other)
   {
      if(this == &other)
      {
         return *this;
      } else
      {
         // Copy object
         *this = other;

         // regenerate Smart pointers
         this->initPOperators();

         return *this;
      }
   }
   
   void AssociatedLegendreOperator::initPOperators()
   {
      // Projection operator from Spectral to phi component of gradient
      // Treat the CSCS case differently (has 2 extraneous grid points)
      if(SimulationConstants::isCSCSGrid())
      {
         this->mpProj2GradPh = SmartZRPOp(new PolynomialOperator<CRPOperator<DynamoComplex> >(DynamoComplex(0.0, this->m()), this->mSin_1Ones, this->mSin_1Poly));
      } else
      {
         this->mpProj2GradPh = SmartZRPOp(new PolynomialOperator<CRPOperator<DynamoComplex> >(DynamoComplex(0.0, this->m()), this->mSin_1Theta, this->poly()));
      }

      // Projection operator from Spectral to theta component of gradient
      this->mpProj2GradTh = SmartPOp(new PolynomialOperator<POperator>(this->diff(1)));

      // Projection operator from Spectral T compoent to phi component
      this->mpProjT2Ph = SmartLPOp(new PolynomialOperator<LPOperator>(this->sll_1(), this->diff(1)));

      // Projection operator from Spectral T compoent to theta component
      // Treat the CSCS case differently (has 2 extraneous grid points)
      if(SimulationConstants::isCSCSGrid())
      {
         this->mpProjT2Th = SmartZLRPOp(new PolynomialOperator<CLRPOperator<DynamoComplex> >(DynamoComplex(0.0, -this->m()), this->sll_1(), this->mSin_1Ones, this->mSin_1Poly));
      } else
      {
         this->mpProjT2Th = SmartZLRPOp(new PolynomialOperator<CLRPOperator<DynamoComplex> >(DynamoComplex(0.0, -this->m()), this->sll_1(), this->mSin_1Theta, this->poly()));
      }

      // Projection operator from Spectral S compoent to phi component
      // Treat the CSCS case differently (has 2 extraneous grid points)
      if(SimulationConstants::isCSCSGrid())
      {
         this->mpProjS2Ph = SmartZLRPOp(new PolynomialOperator<CLRPOperator<DynamoComplex> >(DynamoComplex(0.0, this->m()), this->sll_1(), this->mSin_1Ones, this->mSin_1Poly));
      } else
      {
         this->mpProjS2Ph = SmartZLRPOp(new PolynomialOperator<CLRPOperator<DynamoComplex> >(DynamoComplex(0.0, this->m()), this->sll_1(), this->mSin_1Theta, this->poly()));
      }

      // Projection operator from Spectral S compoent to theta component
      this->mpProjS2Th = SmartLPOp(new PolynomialOperator<LPOperator>(this->sll_1(), this->diff(1)));

      // Projection operator from Spectral to first derivative
      this->mpDProj = SmartPOp(new PolynomialOperator<POperator>(this->diff(1)));

      // Projection operator from Spectral to polynomial
      this->mpProj = SmartPOp(new PolynomialOperator<POperator>(this->poly()));

      // Integration operator from real space phi component to spectral T component
      this->mpIntgPh2T = SmartRPOp(new PolynomialOperator<RPOperator>(this->sll_1(), this->wDiff(1)));

//      // Integration operator from real space theta component to spectral T component
      this->mpIntgTh2T =  SmartZLRPOp(new PolynomialOperator<CLRPOperator<DynamoComplex> >(DynamoComplex(0.0, this->m()), this->mSin_1Theta, this->sll_1(), this->wPoly()));

//      // Integration operator from real space phi component to spectral S component
      this->mpIntgPh2S = SmartZLRPOp(new PolynomialOperator<CLRPOperator<DynamoComplex> >(DynamoComplex(0.0, -this->m()), this->mSin_1Theta, this->sll_1(), this->wPoly()));

      // Integration operator from real space theta component to spectral S component
      this->mpIntgTh2S = SmartRPOp(new PolynomialOperator<RPOperator>(this->sll_1(), this->wDiff(1)));

      // Integration operator from real space to spectral
      this->mpIntg = SmartPOp(new PolynomialOperator<POperator>(this->wPoly()));
   }
}
