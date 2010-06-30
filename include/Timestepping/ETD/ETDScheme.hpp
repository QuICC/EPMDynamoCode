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
          * @param hasL0 Is l=0 mode required?
          */
         ETDScheme(EPMFloat a, EPMFloat b, const BasisType &basis, TimestepParameters &tsteps, SmartTruncation pTrunc, bool hasL0);

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

      private:
   };

   template <typename TSimType, template <typename> class TMethod> ETDScheme<TSimType, TMethod>::ETDScheme(EPMFloat a, EPMFloat b, const typename ETDScheme<TSimType, TMethod>::BasisType &basis, TimestepParameters &tsteps, SmartTruncation pTrunc, bool hasL0)
      : TMethod<TSimType>(a, b, basis, tsteps, pTrunc, hasL0) 
   {
   }

   template <typename TSimType, template <typename> class TMethod> void ETDScheme<TSimType, TMethod>::timestep(typename ETDScheme<TSimType, TMethod>::ScalarType& rVar, typename ETDScheme<TSimType, TMethod>::ScalarType& nTerms)
   {
      // If the timestep has been rejected recover previous timestep values
      if(this->rTSParams().isRejected())
      {
         // Restore value from previous timestep
         this->restoreOld(rVar, nTerms);
      }

      // Do first iteration including recomputation of matrices if required
      if(this->rTSParams().isNextStep())
      {
         // Reset the iteration pointer
         this->resetIterations();

         // Store the variable before timestep to allow rejection of timestep
         this->storeOld(rVar, nTerms);

         // Update timesteping matrices if required
         if(this->rTSParams().hasNewDt())
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
}

#endif // ETDSCHEME_HPP
