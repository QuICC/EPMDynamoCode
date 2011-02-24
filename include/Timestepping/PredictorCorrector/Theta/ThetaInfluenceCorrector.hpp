/** \file ThetaInfluenceCorrector.hpp
 *  \brief Implemenation of the corrector step of the theta method
 */

#ifndef THETAINFLUENCECORRECTOR_HPP
#define THETAINFLUENCECORRECTOR_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "Timestepping/PredictorCorrector/Theta/ThetaCorrector.hpp"
#include "Timestepping/PredictorCorrector/Theta/ThetaLHSTOperatorSet.hpp"
#include "Timestepping/PredictorCorrector/Theta/ThetaRHSTOperatorSet.hpp"
#include "Equations/InfluenceMatrix.hpp"

namespace EPMDynamo {

   /**
    * \brief Implemenation of the corrector step of the theta method
    *
    * \tparam TSimType Type of the simulation
    */
   template <typename TSimType> class ThetaInfluenceCorrector: public ThetaCorrector<TSimType>
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

         /// Typedef for a smart pointer to an influence matrix operator
         typedef EPMSHARED_PTR<InfluenceMatrix<TSimType> > SmartInfluenceMatrix;

         /**
          * @brief Constructor
          *
          * @param pOpLHS Pointer to the LHS operator
          * @param pOpRHS Pointer to the RHS operator
          */
         ThetaInfluenceCorrector(SmartInfluenceMatrix pInfluence, SmartScalarType pPrevious, SmartLHSOperators pOpLHS, SmartRHSOperators pOpRHS);

         /**
          * @brief Destructor
          */
         virtual ~ThetaInfluenceCorrector() {};

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
          * @brief Set the corrector step RHS
          *
          * @param newNTerms New non linear terms
          */
         void setRHS(ScalarType& rNTerms);

         /**
          * @brief Solve the corrector equations
          */
         void solve(ScalarType& rNTerms);

      private:
   };

   template <typename TSimType> ThetaInfluenceCorrector<TSimType>::ThetaInfluenceCorrector(typename ThetaInfluenceCorrector<TSimType>::SmartInfluenceMatrix pInfluence, typename ThetaInfluenceCorrector<TSimType>::SmartScalarType pPrevious, typename ThetaInfluenceCorrector<TSimType>::SmartLHSOperators pOpLHS, typename ThetaInfluenceCorrector<TSimType>::SmartRHSOperators pOpRHS)
      : ThetaCorrector<TSimType>(pPrevious, pOpLHS, pOpRHS), mpInfluence(pInfluence)
   {
   }

   template <typename TSimType> void ThetaInfluenceCorrector<TSimType>::compute(typename ThetaInfluenceCorrector<TSimType>::ScalarType &rVar, typename ThetaInfluenceCorrector<TSimType>::ScalarType &rNTerms)
   {
      // Set the RHS part of corrector step
      this->setRHS(rNTerms);

      // Solve corrector equations
      this->solve(rNTerms);

      // Use corrector solution
      this->useCorrection(rVar, rNTerms);
   }

   template <typename TSimType> void ThetaInfluenceCorrector<TSimType>::setRHS(typename ThetaInfluenceCorrector<TSimType>::ScalarType& rNTerms)
   {
      // Solve influenc matrix equation
      this->mpInfluence->solve(rNTerms);
      
      // Compute RHS part of corrector timestep equation
      ThetaCorrector<TSimType>::setRHS(rNTerms);
   }

   template <typename TSimType> void ThetaInfluenceCorrector<TSimType>::solve(typename ThetaInfluenceCorrector<TSimType>::ScalarType& rNTerms)
   {
      // Get the correction to the unknown variable
      ThetaCorrector<TSimType>::solve(rNTerms);

      // Include kernel influence
      this->mpInfluence->addKernel(rNTerms);
   }
}

#endif // THETAINFLUENCECORRECTOR_HPP
