/** \file ThetaPredictor.hpp
 *  \brief Implemenation of the predictor step of the theta method
 */

#ifndef THETAPREDICTOR_HPP
#define THETAPREDICTOR_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "Timestepping/SchemeIteration.hpp"
#include "Timestepping/ErrorControl.hpp"
#include "Timestepping/PredictorCorrector/Theta/ThetaLHSTOperatorSet.hpp"
#include "Timestepping/PredictorCorrector/Theta/ThetaRHSTOperatorSet.hpp"

namespace EPMDynamo {

   /**
    * \brief Implemenation of the predictor step of the theta method
    *
    * \tparam TSimType Type of the simulation
    */
   template <typename TSimType> class ThetaPredictor: public SchemeIteration<TSimType>
   {
      public:
         /// Typedef from Simulation trait to local scalar type
         typedef typename TSimType::ScalarType    ScalarType;

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
         ThetaPredictor(SmartLHSOperators pOpLHS, SmartRHSOperators pOpRHS);

         /**
          * @brief Destructor
          */
         virtual ~ThetaPredictor() {};

         /**
          * @brief Compute the implemented iteration
          *
          * @param rVar Input/Output variable
          * @param rNTerms New non linear terms
          */
         virtual void compute(ScalarType &rVar, ScalarType &rNTerms);
         
      protected:
         /**
          * @brief Pointer to the LHS operator
          */
         SmartLHSOperators  mpOpLHS;

         /**
          * @brief Pointer to the RHS operator
          */
         SmartRHSOperators  mpOpRHS;

         /**
          * @brief Set the predictor step RHS
          *
          * @param rVar Input/Output variable
          * @param oldVar Previous variable
          * @param nTerms Non linear terms
          */
         void setRHS(ScalarType& rVar, const ScalarType& oldVar, const ScalarType& nTerms);

         /**
          * @brief Solve the predictor equations
          *
          * @param rVar Input/Output variable
          */
         void solve(ScalarType& rVar);

      private:
   };

   template <typename TSimType> ThetaPredictor<TSimType>::ThetaPredictor(SmartLHSOperators pOpLHS, SmartRHSOperators pOpRHS)
      : mpOpLHS(pOpLHS), mpOpRHS(pOpRHS)
   {
   }

   template <typename TSimType> void ThetaPredictor<TSimType>::compute(typename ThetaPredictor<TSimType>::ScalarType &rVar, typename ThetaPredictor<TSimType>::ScalarType &rNTerms)
   {
      // Set the RHS part of predictor step
      this->setRHS(rVar, rVar, rNTerms);

      // Solve predictor equations
      this->solve(rVar);
   }

   template <typename TSimType> void ThetaPredictor<TSimType>::setRHS(typename ThetaPredictor<TSimType>::ScalarType& rVar, const typename ThetaPredictor<TSimType>::ScalarType& oldVar, const typename ThetaPredictor<TSimType>::ScalarType& nTerms)
   {
      // Get number of harmonic degrees
      int nL = this->mpOpRHS->trunc()->local()->spec()->nL();
      const int l0 = rVar.minL();

      // Loop over degrees
      for(int l = l0; l < nL; ++l)
      {
         // Compute RHS part of predictor timestep equation for given degree
         this->mpOpRHS->affineOrders(rVar.rLShell(l), oldVar.lshell(l), nTerms.lshell(l), l);
      }
   }

   template <typename TSimType> void ThetaPredictor<TSimType>::solve(typename ThetaPredictor<TSimType>::ScalarType& rVar)
   {
      // Solve for unknown variable
      this->mpOpLHS->solve(rVar);
   }
}

#endif // THETAPREDICTOR_HPP
