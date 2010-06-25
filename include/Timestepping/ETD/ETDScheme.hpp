/** \file ETDScheme.hpp
 *  \brief Implementation of a general ETD scheme
 */

#ifndef ETDSCHEME_HPP
#define ETDSCHEME_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "Domain/Truncation.hpp"
#include "Simulations/Traits/SimulationTraits.hpp"
#include "Timestepping/TimestepParameters.hpp"

namespace EPMDynamo {

   /**
    * \brief Implementation of a general Predictor/Corrector scheme
    *
    * \tparam TSimType Type of the simulation
    * \tparam TMethod ETD Method to use
    */
   template <typename TSimType, template <typename> class TMethod> class ETDScheme: public TMethod<TSimType>
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
         ETDScheme(EPMFloat a, EPMFloat b, const BasisType &basis, TimestepParameters &tsteps, SmartTruncation pTrunc);

         /**
          * @brief Destructor
          */
         virtual ~ETDScheme() {};

         /**
          * @brief Perform a type step (the actual order/type depends on template)
          *
          * @param rVar Input/Output variable to timestep
          * @param nTerms Non linear terms
          */
         void timestep(ScalarType& rVar, ScalarType& nTerms);
         
      protected:

         /**
          * @brief Compute an iteration of scheme (actuall computation depends on type)
          *
          * @param rVar Input/Output variable
          * @param newNTerms New non linear terms
          */
         void doIteration(ScalarType& rVar, ScalarType& nTerms);

      private:
   };

   template <typename TSimType, template <typename> class TMethod> ETDScheme<TSimType, TMethod>::ETDScheme(EPMFloat a, EPMFloat b, const typename ETDScheme<TSimType, TMethod>::BasisType &basis, TimestepParameters &tsteps, SmartTruncation pTrunc)
      : TMethod<TSimType>(a, b, basis, tsteps, pTrunc) 
   {
   }

   template <typename TSimType, template <typename> class TMethod> void ETDScheme<TSimType, TMethod>::timestep(typename ETDScheme<TSimType, TMethod>::ScalarType& rVar, typename ETDScheme<TSimType, TMethod>::ScalarType& nTerms)
   {
      // If the timestep has been rejected recover previous timestep values
      if(this->mrTStepParams.isRejected())
      {
         // Restore value from previous timestep
         this->restoreOld(rVar, nTerms);
      }

      // Do first iteration including recomputation of matrices if required
      if(this->mrTStepParams.isNextStep())
      {
         // Store the variable before timestep to allow rejection of timestep
         this->storeOld(rVar, nTerms);

         // Update timesteping matrices if required
         if(this->mrTStepParams.hasNewDt())
         {
            this->updateTimeMatrices();
         }

         // Do iteration (the actual work done will depend on precise method)
         this->doIteration(rVar, nTerms);
      } else
      {
         // Do next iteration (the actual work done will depend on precise method)
         this->doIteration(rVar, nTerms);
      }
   }

   template <typename TSimType, template <typename> class TMethod> void ETDScheme<TSimType, TMethod>::doIteration(typename ETDScheme<TSimType, TMethod>::ScalarType& rVar, typename ETDScheme<TSimType, TMethod>::ScalarType& nTerms)
   {
      // Check for intermediate computations
      if(this->hasIntermediate())
      {
         // Do intermediate step
         this->computeIntermediate(rVar, nTerms);
      } else
      {
         // Do finale timestep
         this->doStep(rVar, nTerms);
      }
   }

}

#endif // ETDSCHEME_HPP
