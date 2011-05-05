/** \file ThetaMethod.hpp
 *  \brief Implementation of the predictor-corrector theta method (without influence matrix)
 */

#ifndef THETAMETHOD_HPP
#define THETAMETHOD_HPP

// Configuration includes
//
#include "Config/SmartPointer.h"

// System includes
//

// External includes
//

// Project includes
//
#include "Domain/Truncation.hpp"
#include "Timestepping/TimestepParameters.hpp"
#include "Timestepping/IterativeSchemeBase.hpp"
#include "Timestepping/PredictorCorrector/Theta/ThetaLHSTOperatorSet.hpp"
#include "Timestepping/PredictorCorrector/Theta/ThetaRHSTOperatorSet.hpp"
#include "Timestepping/PredictorCorrector/Theta/ThetaPredictor.hpp"
#include "Timestepping/PredictorCorrector/Theta/ThetaCorrector.hpp"
#include "Config/SimulationConfig.hpp"

namespace EPMDynamo {

   /**
    * \brief Implementation of the predictor-corrector theta method (without influence matrix)
    *
    * \tparam TSimType Type of the simulation
    */
   template <typename TSimType> class ThetaMethod: public IterativeSchemeBase<TSimType>
   {
      public:
         /// Typedef from Simulation trait to local radial basis type
         typedef typename TSimType::RadialBasisType    BasisType;

         /// Typedef from Simulation trait to local scalar type
         typedef typename TSimType::ScalarType    ScalarType;

         /// Typedef for a smart pointer to LHS operator set
         typedef EPMSHARED_PTR<ThetaLHSTOperatorSet<TSimType> > SmartLHSOperator;

         /// Typedef for a smart pointer to RHS operator set
         typedef EPMSHARED_PTR<ThetaRHSTOperatorSet<TSimType> > SmartRHSOperator;

         /**
          * @brief Constructor
          *
          * @param a Coefficient \f$a\f$ of timestep scheme
          * @param b Coefficient \f$b\f$ of timestep scheme
          * @param basis Reference to the basis used for the operators
          * @param tsteps Timestep parameters
          * @param pTrunc Truncation information
          * @param hasL0 Is l=0 mode required?
          */
         ThetaMethod(EPMFloat a, EPMFloat b, const BasisType &basis, TimestepParameters &tsteps, SmartTruncation pTrunc, bool hasL0);

         /**
          * @brief Destructor
          */
         virtual ~ThetaMethod() {};

         /**
          * @brief Add a boundary condition
          *
          * @param pBC Boundary condition
          */
         virtual void addBC(SmartBC pBC);

         /**
          * @brief Initialise the method
          */
         void init();
         
      protected:
         /**
          * @brief Timestepping matrices for the LHS part of equation
          *
          * These matrices are the LU decomposition of the timestep operators
          */
         SmartLHSOperator  mpLHS;

         /**
          * @brief Timestepping matrices for the RHS part of equation
          */
         SmartRHSOperator  mpRHS;

         /**
          * @brief Update the timestep matrices after a timestep change
          */
         void updateTimeMatrices();

         /**
          * @brief Initialise the operators
          */
         void initOperators();

         /**
          * @brief Compute next scheme iteration
          *
          * @param rVar Input/Output variable
          * @param rNTerms New non linear terms
          */
         void doIteration(ScalarType& rVar, ScalarType& rNTerms);

      private:

         /**
          * @brief Initialise the storage for the method
          */
         void initStorage();

         /**
          * @brief Update values of previous iteration NTerms
          */
         void updatePrevious(const ScalarType& previous);
   };

   template <typename TSimType> ThetaMethod<TSimType>::ThetaMethod(EPMFloat a, EPMFloat b, const typename ThetaMethod<TSimType>::BasisType &basis, TimestepParameters &tsteps, SmartTruncation pTrunc, bool hasL0)
      : IterativeSchemeBase<TSimType>(a, b, basis, tsteps, pTrunc, hasL0), mpLHS(new ThetaLHSTOperatorSet<TSimType>(a, b, basis, pTrunc, hasL0)), mpRHS(new ThetaRHSTOperatorSet<TSimType>(a, b, basis, pTrunc, hasL0))
   {
   }

   template <typename TSimType> void ThetaMethod<TSimType>::addBC(SmartBC pBC)
   {
      // Add boundary condition to LHS of theta method
      this->mpLHS->addBC(pBC);

      // Add boundary condition to RHS of theta method
      this->mpRHS->addBC(pBC);
   }

   template <typename TSimType> void ThetaMethod<TSimType>::init()
   {
      // initialise pointers
      this->initStorage();

      // initialise the operators
      this->initOperators();
   }

   template <typename TSimType> void ThetaMethod<TSimType>::initOperators()
   {
      // Create LHS operators and LU factorise them
      this->mpLHS->initOperators();

      // Create RHS operators
      this->mpRHS->initOperators();
   }

   template <typename TSimType> void ThetaMethod<TSimType>::updateTimeMatrices()
   {
      // Updated the LHS operator due to new timestep
      this->mpLHS->update(this->rTSParams().dt());

      // Updated the RHS operator due to new timestep
      this->mpRHS->update(this->rTSParams().dt());
   }

   template <typename TSimType> void ThetaMethod<TSimType>::initStorage()
   {
      //
      // Create intermediate storage
      EPMSHARED_PTR<ScalarType > pPrevious(new ScalarType (this->mpLHS->trunc(), this->mpLHS->hasL0()));
      this->mTmpScalars.push_back(pPrevious);

      // Create the predictor computation step
      EPMSHARED_PTR<ThetaPredictor<TSimType> > pItP(new ThetaPredictor<TSimType> (this->mpLHS, this->mpRHS));

      // Create the corrector computation step
      EPMSHARED_PTR<ThetaCorrector<TSimType> > pItC(new ThetaCorrector<TSimType> (pPrevious, this->mpLHS, this->mpRHS));

      //
      // Add the required computation steps
      
      // Add the predictor computation
      this->mSteps.push_back(pItP);
      // Add a corrector step
      this->mSteps.push_back(pItC);
   }

   template <typename TSimType> void ThetaMethod<TSimType>::doIteration(typename ThetaMethod<TSimType>::ScalarType& rVar, typename ThetaMethod<TSimType>::ScalarType& rNTerms)
   {
      // Do normal iteration
      IterativeSchemeBase<TSimType>::doIteration(rVar, rNTerms);

      // Update values stored from previous iteration
      this->updatePrevious(rNTerms);
   }

   template <typename TSimType> inline void ThetaMethod<TSimType>::updatePrevious(const typename ThetaMethod::ScalarType& previous)
   {
      // Get number of harmonic degrees
      int nL = this->mTmpScalars.at(0)->nL();
      // Get minimu harmonic degrees
      const int l0 = this->mTmpScalars.at(0)->minL();

      // Loop over degrees
      for(int l = l0; l < nL; ++l)
      {
         // Store previous values
         this->mTmpScalars.at(0)->rLShell(l) = previous.lshell(l);
      }
   }

}

#endif // THETAMETHOD_HPP
