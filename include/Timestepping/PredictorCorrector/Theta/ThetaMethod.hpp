/** \file ThetaMethod.hpp
 *  \brief Implementation of the predictor-corrector theta method (without influence matrix)
 */

#ifndef THETAMETHOD_HPP
#define THETAMETHOD_HPP

// Configuration includes
//
#include "Config/SmartPointer.h"
#include "Config/SimulationConfig.hpp"

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

namespace EPMDynamo {

   /**
    * \brief Implementation of the predictor-corrector theta method (without influence matrix)
    */
   class ThetaMethod: public IterativeSchemeBase
   {
      public:
         /// Typedef from Simulation trait to local radial basis type
         typedef SimulationConfig::NumericalScheme::RadialBasisType    BasisType;

         /// Typedef from Simulation trait to local scalar type
         typedef SimulationConfig::NumericalScheme::ScalarType    ScalarType;

         /// Typedef for a smart pointer to LHS operator set
         typedef EPMSHARED_PTR<ThetaLHSTOperatorSet > SmartLHSOperator;

         /// Typedef for a smart pointer to RHS operator set
         typedef EPMSHARED_PTR<ThetaRHSTOperatorSet > SmartRHSOperator;

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

   inline void ThetaMethod::updatePrevious(const ThetaMethod::ScalarType& previous)
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
