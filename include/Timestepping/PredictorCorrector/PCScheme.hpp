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
#include "Simulations/Traits/SimulationTraits.hpp"
#include "Timestepping/TimestepParameters.hpp"

namespace EPMDynamo {

   /**
    * \brief Implementation of a general Predictor/Corrector scheme
    *
    * \tparam TSimType Type of the simulation
    * \tparam TMethod PC Method used
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
          * @param hasL0 Is l=0 mode required?
          */
         PCScheme(EPMFloat a, EPMFloat b, const BasisType &basis, TimestepParameters &tsteps, SmartTruncation pTrunc, bool hasL0);

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
         
      protected:
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

   template <typename TSimType, template <typename> class TMethod> PCScheme<TSimType, TMethod>::PCScheme(EPMFloat a, EPMFloat b, const typename PCScheme<TSimType, TMethod>::BasisType &basis, TimestepParameters &tsteps, SmartTruncation pTrunc, bool hasL0)
      : TMethod<TSimType>(a, b, basis, tsteps, pTrunc, hasL0)
   {
   }

   template <typename TSimType, template <typename> class TMethod> void PCScheme<TSimType, TMethod>::timestep(typename PCScheme<TSimType, TMethod>::ScalarType& rVar, typename PCScheme<TSimType, TMethod>::ScalarType& nTerms)
   {
      // If the timestep has been rejected recover previous timestep values
      if(this->rTSParams().isRejected())
      {
         // Restore value from previous timestep
         this->restoreOld(rVar, nTerms);
      }

      // Do predictor step if previous corrector solution converged
      if(this->rTSParams().isNextStep())
      {
         // Store the variable before timestep to allow rejection of timestep
         this->storeOld(rVar, nTerms);

         // Update timesteping matrices if required
         if(this->rTSParams().hasNewDt())
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
