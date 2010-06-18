/** \file PCScheme.hpp
 *  \brief Implementation of a general Predictor/Corrector scheme
 */

#ifndef PCSCHEME_HPP
#define PCSCHEME_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "Domain/Truncation.hpp"
#include "BoundaryConditions/BoundaryCondition.hpp"
#include "Simulations/SimulationTraits.hpp"
#include "Timestepping/TimestepParameters.hpp"

namespace EPMDynamo {

   /**
    * \brief Implementation of a general Predictor/Corrector scheme
    *
    * \tparam TSimType Type of the simulation
    * \tparam TMethod PC Method used
    *
    * \bug Needs big cleaning and restructuration
    */
   template <typename TSimType, template <typename> class TMethod> class PCScheme: public TMethod<TSimType>
   {
      public:
         /// Typedef from Simulation trait to local radial basis type
         typedef typename TSimType::RadialBasisType    BasisType;

         /// Typedef from Simulation trait to local scalar type
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
         PCScheme(DynamoFloat a, DynamoFloat b, const BasisType &basis, TimestepParameters &tsteps, SmartTruncation pTrunc);

         /**
          * @brief Simple empty destructor
          */
         virtual ~PCScheme() {};

         /**
          * @brief Perform a type step (either predictor or corrector)
          *
          * @param rVar Input/Output variable to timestep
          * @param nTerms Non linear terms
          */
         void timestep(ScalarType& rVar, ScalarType& nTerms);

         /**
          * @brief Set the theta parameter of the scheme
          *
          * @param theta Set implicitness
          *
          * \bug Implicitness implementation and handling has to be cleaned!!
          */
         void setTheta(DynamoFloat theta);
         
      protected:
         /**
          * @brief Values of the variable of previous step
          *
          * \bug TO BE REMOVED!!
          */
         ScalarType  mTimeDiff;

         /**
          * @brief Compute the temporal truncation errors
          *
          * @param var Input variable
          *
          * \bug TO BE REMOVED!!
          */
         void computeTruncationError(const ScalarType& var);

         /**
          * @brief Do predictor timestep
          *
          * @param rVar Input/output variable
          * @param oldVar Previous values
          * @param nTerms Non linear terms
          */
         void doPredictorStep(ScalarType& rVar, const ScalarType& oldVar, ScalarType& nTerms);

         /**
          * @brief Do corrector timestep
          *
          * @param rVar Input/Output variable
          * @param newNTerms New non linear terms
          */
         void doCorrectorStep(ScalarType& rVar, ScalarType& newNTerms);

      private:
   };

   template <typename TSimType, template <typename> class TMethod> PCScheme<TSimType, TMethod>::PCScheme(DynamoFloat a, DynamoFloat b, const typename PCScheme<TSimType, TMethod>::BasisType &basis, TimestepParameters &tsteps, SmartTruncation pTrunc)
      : TMethod<TSimType>(a, b, basis, tsteps, pTrunc), mTimeDiff(pTrunc) 
   {
   }

   template <typename TSimType, template <typename> class TMethod> void PCScheme<TSimType, TMethod>::setTheta(DynamoFloat theta)
   {
      // Set theta parameter in LHS operators
      this->mLHS.setTheta(theta);

      // Set theta parameter in RHS operators
      this->mRHS.setTheta(theta);
   }

   template <typename TSimType, template <typename> class TMethod> void PCScheme<TSimType, TMethod>::computeTruncationError(const typename PCScheme<TSimType, TMethod>::ScalarType& var)
   {
      // Get number of harmonic degrees
      int nL = this->oldNTerms().nL();
      // Get minimu harmonic degrees
      const int l0 = this->oldNTerms().minL();

      DynamoFloat err = 0.0;
      DynamoFloat err2 = 0.0;
      DynamoFloat dt = this->mrTStepParams.dt();

      // Loop over degrees
      for(int l = l0; l < nL; ++l)
      {
         // Compute norm of the second time derivative approximation
         err2 += 0.5*dt*((this->mTimeDiff.lshell(l) - (this->oldVar().lshell(l) - var.lshell(l))/dt).norm());

         // First time derivative approximation
         this->mTimeDiff.rLShell(l) = (this->oldVar().lshell(l) - var.lshell(l))/dt;

         // Compute first derivative error approximation
         err += dt*this->mTimeDiff.lshell(l).norm();
      }

      // Store first derivative error
      this->mrTStepParams.updateDtError(1, err);

      // Store second derivative error
      this->mrTStepParams.updateDtError(2, err2);
   }

   template <typename TSimType, template <typename> class TMethod> void PCScheme<TSimType, TMethod>::timestep(typename PCScheme<TSimType, TMethod>::ScalarType& rVar, typename PCScheme<TSimType, TMethod>::ScalarType& nTerms)
   {
      // If the timestep has been rejected recover previous timestep values
      if(this->mrTStepParams.isRejected())
      {
         // Restore value from previous timestep
         this->restoreOld(rVar, nTerms);
      }

      // Do predictor step if previous corrector solution converged
      if(this->mrTStepParams.isNextStep())
      {
// This is only here temporary until error output is not required anymore
this->mrTStepParams.resetDtError(1);
this->mrTStepParams.resetDtError(2);
this->mrTStepParams.resetError();

         // Store the variable before timestep to allow rejection of timestep
         this->storeOld(rVar, nTerms);

         // Update timesteping matrices if required
         if(this->mrTStepParams.hasNewDt())
         {
            this->updateTimeMatrices();
         }

         // Do predictor step
         this->doPredictorStep(rVar, rVar, nTerms);

      // Do corrector step if previous corrector iteration solution is not converged yet
      } else
      {
         // Do corrector step
         this->doCorrectorStep(rVar, nTerms);

         // Compute the temporal truncation error
         this->computeTruncationError(rVar);
      }

      // Update values stored from previous iteration
      this->updatePrevious(nTerms);
   }

   template <typename TSimType, template <typename> class TMethod> void PCScheme<TSimType, TMethod>::doPredictorStep(typename PCScheme<TSimType, TMethod>::ScalarType& rVar, const typename PCScheme<TSimType, TMethod>::ScalarType& oldVar, typename PCScheme<TSimType, TMethod>::ScalarType& nTerms)
   {
      // Set the RHS part of predictor step
      this->setPredictorRHS(rVar, oldVar, nTerms);

      // Solve predictor equations
      this->solvePredictor(rVar);
   }

   template <typename TSimType, template <typename> class TMethod> void PCScheme<TSimType, TMethod>::doCorrectorStep(typename PCScheme<TSimType, TMethod>::ScalarType& rVar, typename PCScheme<TSimType, TMethod>::ScalarType& newNTerms)
   {
      // Set the RHS part of corrector step
      this->setCorrectorRHS(newNTerms);

      // Solve corrector equations
      this->solveCorrector();

      // Use corrector solution
      this->useCorrection(rVar);
   }
}

#endif // PCSCHEME_HPP
