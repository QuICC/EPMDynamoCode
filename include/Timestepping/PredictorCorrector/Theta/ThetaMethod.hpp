/** \file ThetaMethod.hpp
 *  \brief Implementation of the theta method (without influence matrix)
 */

#ifndef THETAMETHOD_HPP
#define THETAMETHOD_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "Domain/Truncation.hpp"
#include "BoundaryConditions/BoundaryCondition.hpp"
#include "Simulations/Traits/SimulationTraits.hpp"
#include "Timestepping/TimestepSchemeBase.hpp"
#include "Timestepping/PredictorCorrector/PCErrorControl.hpp"
#include "Timestepping/PredictorCorrector/Theta/ThetaLHSTOperatorSet.hpp"
#include "Timestepping/PredictorCorrector/Theta/ThetaRHSTOperatorSet.hpp"
#include "Timestepping/TimestepParameters.hpp"

namespace EPMDynamo {

   /**
    * \brief Implementation of the theta method (without influence matrix)
    *
    * \tparam TSimType Type of the simulation
    *
    * \bug Needs to be reviewd and implemented in a clean way
    */
   template <typename TSimType> class ThetaMethod: public TimestepSchemeBase<TSimType>
   {
      public:
         /// Typedef from Simulation trait to local truncation type
         typedef typename TSimType::RadialBasisType    BasisType;

         /// Typedef from Simulation trait to local truncation type
         typedef typename TSimType::ScalarType    ScalarType;

         /**
          * @brief Constructor
          *
          * @param a Coefficient \f$a\f$ of timestep scheme
          * @param b Coefficient \f$b\f$ of timestep scheme
          * @param basis Reference to the basis used for the operators
          * @param tsteps Timestep parameters
          * @param pTrunc Truncation information
          */
         ThetaMethod(EPMFloat a, EPMFloat b, const BasisType &basis, TimestepParameters &tsteps, SmartTruncation pTrunc);

         /**
          * @brief Simple empty destructor
          */
         virtual ~ThetaMethod() {};

         /**
          * @brief Add a boundary condition
          *
          * @param pBC Boundary condition
          */
         void addBC(SmartBC pBC);

         /**
          * @brief Initialise the ThetaMethod
          */
         void init();
         
      protected:
         /**
          * @brief Timestepping matrices for the LHS part of equation
          *
          * These matrices are the LU decomposition of the timestep operators
          */
         ThetaLHSTOperatorSet<TSimType>  mLHS;

         /**
          * @brief Timestepping matrices for the RHS part of equation
          */
         ThetaRHSTOperatorSet<TSimType>  mRHS;

         /**
          * @brief Storage for the values of the previous iteration
          */
         ScalarType  mPreviousNTerms;

         /**
          * @brief Update values of previous iteration
          *
          * @param previous Previous variable
          */
         void updatePrevious(const ScalarType& previous);

         /**
          * @brief Update the timestepping operators after change of timestep length
          */
         void updateTimeMatrices();

         /**
          * @brief Set the predictor step RHS
          *
          * @param rVar Input/Output variable
          * @param oldVar Previous variable
          * @param nTerms Non linear terms
          */
         void setPredictorRHS(ScalarType& rVar, const ScalarType& oldVar, const ScalarType& nTerms);

         /**
          * @brief Solve the predictor equations
          *
          * @param rVar Input/Output variable
          */
         void solvePredictor(ScalarType& rVar);

         /**
          * @brief Set the predictor step RHS
          *
          * @param newNTerms New non linear terms
          */
         void setCorrectorRHS(const ScalarType& newNTerms);

         /**
          * @brief Solve the corrector equations
          */
         void solveCorrector();

         /**
          * @brief Use the corrector solution
          *
          * @param rVar Input/Output variable
          *
          * \bug Corrector/Correction norm computatin have to be tested and implemented in a clean way
          */
         void useCorrection(ScalarType& rVar);

         /**
          * @brief Add correction from corrector step
          *
          * @param rVar Input/Output variable
          * @param corr Corrector correction
          */
         void addCorrection(ScalarType& rVar, const ScalarType& corr) const;

      private:
   };

   template <typename TSimType> ThetaMethod<TSimType>::ThetaMethod(EPMFloat a, EPMFloat b, const typename ThetaMethod<TSimType>::BasisType &basis, TimestepParameters &tsteps, SmartTruncation pTrunc)
      : TimestepSchemeBase<TSimType>(tsteps, pTrunc), mLHS(a, b, basis, pTrunc), mRHS(a, b, basis, pTrunc), mPreviousNTerms(pTrunc)
   {
   }

   template <typename TSimType> void ThetaMethod<TSimType>::addBC(SmartBC pBC)
   {
      // Add boundary condition to LHS ThetaMethod
      this->mLHS.addBC(pBC);

      // Add boundary condition to RHS ThetaMethod
      this->mRHS.addBC(pBC);
   }

   template <typename TSimType> void ThetaMethod<TSimType>::init()
   {
      // Set index controler for LHS and create operators (force creation of all homogeneous operators because of corrector step)
      this->mLHS.initOperators();

      // Set index controler for RHS
      this->mRHS.initOperators();
   }

   template <typename TSimType> void ThetaMethod<TSimType>::updateTimeMatrices()
   {
      // Updated the LHS operator due to new timestep
      this->mLHS.update(this->rTSParams().dt());

      // Updated the RHS operator due to new timestep
      this->mRHS.update(this->rTSParams().dt());
   }

   template <typename TSimType> inline void ThetaMethod<TSimType>::updatePrevious(const typename ThetaMethod::ScalarType& previous)
   {
      // Get number of harmonic degrees
      int nL = this->mPreviousNTerms.nL();
      // Get minimu harmonic degrees
      const int l0 = this->mPreviousNTerms.minL();

      // Loop over degrees
      for(int l = l0; l < nL; ++l)
      {
         // Store previous values
         this->mPreviousNTerms.rLShell(l) = previous.lshell(l);
      }
   }

   template <typename TSimType> void ThetaMethod<TSimType>::setPredictorRHS(typename ThetaMethod<TSimType>::ScalarType& rVar, const typename ThetaMethod<TSimType>::ScalarType& oldVar, const typename ThetaMethod<TSimType>::ScalarType& nTerms)
   {
      // Get number of harmonic degrees
      int nL = this->mRHS.trunc()->local()->spec()->nL();
      const int l0 = rVar.minL();

      // Loop over degrees
      for(int l = l0; l < nL; ++l)
      {
         // Compute RHS part of predictor timestep equation for given degree
         this->mRHS.affineOrders(rVar.rLShell(l), oldVar.lshell(l), nTerms.lshell(l), l);
      }
   }

   template <typename TSimType> void ThetaMethod<TSimType>::solvePredictor(typename ThetaMethod<TSimType>::ScalarType& rVar)
   {
      // Solve for unknown variable
      this->mLHS.solve(rVar);
   }

   template <typename TSimType> void ThetaMethod<TSimType>::setCorrectorRHS(const typename ThetaMethod<TSimType>::ScalarType& newNTerms)
   {
      // Get number of harmonic degrees
      int nL = this->mRHS.trunc()->local()->spec()->nL();
      const int l0 = this->mPreviousNTerms.minL();

      // Loop over degrees
      for(int l = l0; l < nL; ++l)
      {
         this->mPreviousNTerms.rLShell(l) = this->mRHS.theta()*(newNTerms.lshell(l) - this->mPreviousNTerms.lshell(l));
      }
   }

   template <typename TSimType> void ThetaMethod<TSimType>::solveCorrector()
   {
      // Get the correction to the unknown variable
      this->mLHS.homogeneousBCSolve(this->mPreviousNTerms);
   }

   template <typename TSimType> void ThetaMethod<TSimType>::useCorrection(typename ThetaMethod<TSimType>::ScalarType& rVar)
   {
      // Update simulation wide value of corrector norm
      this->rTSParams().updateError(PCErrorControl<TSimType>::errorNorm(this->mPreviousNTerms, this->rTSParams().error()));

      // Add correction to unknown
      this->addCorrection(rVar, this->mPreviousNTerms);
   }

   template <typename TSimType> void ThetaMethod<TSimType>::addCorrection(typename ThetaMethod<TSimType>::ScalarType& rVar, const typename ThetaMethod<TSimType>::ScalarType& corr) const
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

#endif // THETAMETHOD_HPP
