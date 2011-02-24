/** \file ThetaInfluencePredictor.hpp
 *  \brief Implemenation of the predictor step of the theta method
 */

#ifndef THETAINFLUENCePREDICTOR_HPP
#define THETAINFLUENCEPREDICTOR_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "Timestepping/PredictorCorrector/Theta/ThetaPredictor.hpp"
#include "Timestepping/PredictorCorrector/Theta/ThetaLHSTOperatorSet.hpp"
#include "Timestepping/PredictorCorrector/Theta/ThetaRHSTOperatorSet.hpp"
#include "Equations/InfluenceMatrix.hpp"

namespace EPMDynamo {

   /**
    * \brief Implemenation of the predictor step of the theta method
    *
    * \tparam TSimType Type of the simulation
    */
   template <typename TSimType> class ThetaInfluencePredictor: public ThetaPredictor<TSimType>
   {
      public:
         /// Typedef from Simulation trait to local scalar type
         typedef typename TSimType::ScalarType    ScalarType;

         /// Typedef for a smart pointer to LHS operator set
         typedef EPMSHARED_PTR<ThetaLHSTOperatorSet<TSimType> > SmartLHSOperators;

         /// Typedef for a smart pointer to RHS operator set
         typedef EPMSHARED_PTR<ThetaRHSTOperatorSet<TSimType> > SmartRHSOperators;

         /// Typedef for a smart pointer to an influence matrix operator
         typedef EPMSHARED_PTR<InfluenceMatrix<TSimType> > SmartInfluenceMatrix;

         /**
          * @brief Constructor
          *
          * @param pOpLHS Pointer to the LHS operator
          * @param pOpRHS Pointer to the RHS operator
          */
         ThetaInfluencePredictor(SmartInfluenceMatrix pInfluence, SmartLHSOperators pOpLHS, SmartRHSOperators pOpRHS);

         /**
          * @brief Destructor
          */
         virtual ~ThetaInfluencePredictor() {};

         /**
          * @brief Compute the implemented iteration
          *
          * @param rVar Input/Output variable
          * @param rNTerms New non linear terms
          */
         virtual void compute(ScalarType &rVar, ScalarType &rNTerms);
         
      protected:
         /**
          * @brief Pointer to the influence matrix operator
          */
         SmartInfluenceMatrix  mpInfluence;

         /**
          * @brief Set the predictor step RHS
          *
          * @param rVar Input/Output variable
          * @param oldVar Previous variable
          * @param nTerms Non linear terms
          */
         void setRHS(ScalarType& rVar, const ScalarType& oldVar, ScalarType& rNTerms);

         /**
          * @brief Solve the predictor equations
          *
          * @param rVar Input/Output variable
          */
         void solve(ScalarType& rVar);

      private:
   };

   template <typename TSimType> ThetaInfluencePredictor<TSimType>::ThetaInfluencePredictor(SmartInfluenceMatrix pInfluence, SmartLHSOperators pOpLHS, SmartRHSOperators pOpRHS)
      : ThetaPredictor<TSimType>(pOpLHS, pOpRHS), mpInfluence(pInfluence)
   {
   }

   template <typename TSimType> void ThetaInfluencePredictor<TSimType>::compute(typename ThetaInfluencePredictor<TSimType>::ScalarType &rVar, typename ThetaInfluencePredictor<TSimType>::ScalarType &rNTerms)
   {
      // Set the RHS part of predictor step
      this->setRHS(rVar, rVar, rNTerms);

      // Solve predictor equations
      this->solve(rVar);
   }

   template <typename TSimType> void ThetaInfluencePredictor<TSimType>::setRHS(typename ThetaInfluencePredictor<TSimType>::ScalarType& rVar, const typename ThetaInfluencePredictor<TSimType>::ScalarType& oldVar, typename ThetaInfluencePredictor<TSimType>::ScalarType& rNTerms)
   {
      // Solve influence matrix equation
      this->mpInfluence->solve(rNTerms);

      // Compute RHS part of predictor timestep equation
      ThetaPredictor<TSimType>::setRHS(rVar, oldVar, rNTerms);
   }

   template <typename TSimType> void ThetaInfluencePredictor<TSimType>::solve(typename ThetaInfluencePredictor<TSimType>::ScalarType& rVar)
   {
      // Solve for unknown variable
      ThetaPredictor<TSimType>::solve(rVar);

      // Include kernel influence
      this->mpInfluence->addKernel(rVar);
   }
}

#endif // THETAINFLUENCEPREDICTOR_HPP
