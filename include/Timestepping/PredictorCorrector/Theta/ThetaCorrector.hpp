/** \file ThetaCorrector.hpp
 *  \brief Implemenation of the corrector step of the theta method
 */

#ifndef THETACORRECTOR_HPP
#define THETACORRECTOR_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "Timestepping/SchemeIteration.hpp"
#include "Timestepping/PredictorCorrector/Theta/ThetaLHSTOperatorSet.hpp"
#include "Timestepping/PredictorCorrector/Theta/ThetaRHSTOperatorSet.hpp"

namespace EPMDynamo {

   /**
    * \brief Implemenation of the corrector step of the theta method
    *
    * \tparam TSimType Type of the simulation
    */
   template <typename TSimType> class ThetaCorrector: public SchemeIteration<TSimType>
   {
      public:
         /// Typedef from Simulation trait to local scalar type
         typedef typename TSimType::ScalarType    ScalarType;

         /// typedef for a pointer to a scalar type
          typedef EPMSHARED_PTR<ScalarType> SmartScalarType;

         /// Typedef for a smart pointer to LHS operator set
         typedef EPMSHARED_PTR<ThetaLHSTOperatorSet<TSimType> > SmartLHSOperators;

         /// Typedef for a smart pointer to RHS operator set
         typedef EPMSHARED_PTR<ThetaRHSTOperatorSet<TSimType> > SmartRHSOperators;

         /**
          * @brief Constructor
          *
          * @param pOpLHS Pointer to the LHS operator
          * @param pOpRHS Pointer to the RHS operator
          */
         ThetaCorrector(SmartScalarType pPrevious, SmartLHSOperators pOpLHS, SmartRHSOperators pOpRHS);

         /**
          * @brief Destructor
          */
         virtual ~ThetaCorrector() {};

         /**
          * @brief Compute the implemented iteration
          *
          * @param rVar Input/Output variable
          * @param rNTerms New non linear terms
          */
         virtual void compute(ScalarType &rVar, ScalarType &rNTerms);

         /**
          * @brief Does this iteration provide timestep error information
          */
         bool providesError() const;

         /**
          * @brief Keep the system in the corrector iteration
          */
         bool doNextIteration() const;
         
      protected:
         /**
          * @brief Pointer to the previous NTerms scalar value
          */
         SmartScalarType   mpPreviousNTerms;

         /**
          * @brief Pointer to the LHS operator
          */
         SmartLHSOperators  mpOpLHS;

         /**
          * @brief Pointer to the RHS operator
          */
         SmartRHSOperators  mpOpRHS;

         /**
          * @brief Set the corrector step RHS
          *
          * @param newNTerms New non linear terms
          */
         void setRHS(ScalarType& rNTerms);

         /**
          * @brief Solve the corrector equations
          */
         void solve(ScalarType& rNTerms);

         /**
          * @brief Use the corrector solution
          *
          * @param rVar Input/Output variable
          */
         void useCorrection(ScalarType& rVar, const ScalarType& nTerms);

         /**
          * @brief Add correction from corrector step
          *
          * @param rVar Input/Output variable
          * @param corr Corrector correction
          */
         void addCorrection(ScalarType& rVar, const ScalarType& corr) const;

      private:
   };

   template <typename TSimType> ThetaCorrector<TSimType>::ThetaCorrector(typename ThetaCorrector<TSimType>::SmartScalarType pPrevious, typename ThetaCorrector<TSimType>::SmartLHSOperators pOpLHS, typename ThetaCorrector<TSimType>::SmartRHSOperators pOpRHS)
      : mpPreviousNTerms(pPrevious), mpOpLHS(pOpLHS), mpOpRHS(pOpRHS)
   {
   }

   template <typename TSimType> inline bool ThetaCorrector<TSimType>::providesError() const
   {
      return true;
   }

   template <typename TSimType> inline bool ThetaCorrector<TSimType>::doNextIteration() const
   {
      return false;
   }

   template <typename TSimType> void ThetaCorrector<TSimType>::compute(typename ThetaCorrector<TSimType>::ScalarType &rVar, typename ThetaCorrector<TSimType>::ScalarType &rNTerms)
   {
      // Set the RHS part of corrector step
      this->setRHS(rNTerms);

      // Solve corrector equations
      this->solve(rNTerms);

      // Use corrector solution
      this->useCorrection(rVar, rNTerms);
   }

   template <typename TSimType> void ThetaCorrector<TSimType>::setRHS(typename ThetaCorrector<TSimType>::ScalarType& rNTerms)
   {
      // Get number of harmonic degrees
      int nL = this->mpOpRHS->trunc()->local()->spec()->nL();
      const int l0 = this->mpPreviousNTerms->minL();

      // Loop over degrees
      for(int l = l0; l < nL; ++l)
      {
         rNTerms.rLShell(l) = ThetaTraits<TSimType>::theta*(rNTerms.lshell(l) - this->mpPreviousNTerms->lshell(l));
      }
   }

   template <typename TSimType> void ThetaCorrector<TSimType>::solve(typename ThetaCorrector<TSimType>::ScalarType& rNTerms)
   {
      // Get the correction to the unknown variable
      this->mpOpLHS->solveZero(rNTerms);
   }

   template <typename TSimType> void ThetaCorrector<TSimType>::useCorrection(typename ThetaCorrector<TSimType>::ScalarType& rVar, const typename ThetaCorrector<TSimType>::ScalarType& nTerms)
   {
      // Add correction to unknown
      this->addCorrection(rVar, nTerms);
   }

   template <typename TSimType> void ThetaCorrector<TSimType>::addCorrection(typename ThetaCorrector<TSimType>::ScalarType& rVar, const typename ThetaCorrector<TSimType>::ScalarType& corr) const
   {
      // Get number of harmonic degrees and minimal index
      int nL = rVar.nL();
      int l0 = rVar.minL();

      // loop over degrees
      for(int l = l0; l < nL; ++l)
      {
         // Add corrector correction to solution
         rVar.rLShell(l) += corr.lshell(l);
      }
   }
}

#endif // THETACORRECTOR_HPP
