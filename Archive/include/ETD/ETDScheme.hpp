/** \file ETDScheme.hpp
 */

#ifndef ETDSCHEME_HPP
#define ETDSCHEME_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "Domain/SpectralTruncation.hpp"
#include "Simulations/SimulationTraits.hpp"
#include "Timestepping/TimestepParameters.hpp"

namespace EPMDynamo {

   template <typename TSim, template <typename> class TMethod> class ETDScheme: public TMethod<TSim>
   {
      public:
         /// Typedef from Simulation trait to local radial basis type
         typedef typename TSim::RadialBasisType    BasisType;

         /// Typedef from Simulation trait to local scalar type
         typedef typename TSim::ScalarType    ScalarType;

         /**
          * @brief Constructor
          *
          * @param a Coefficient \f$a\f$ of timestep scheme
          * @param b Coefficient \f$b\f$ of timestep scheme
          * @param basis Reference to the basis used for the operators
          */
         ETDScheme(DynamoFloat a, DynamoFloat b, const BasisType &basis, TimestepParameters &tsteps, SmartSTrunc pTrunc);

         /**
          * @brief Destructor
          */
         virtual ~ETDScheme();

         /**
          * @brief Perform a type step (the actual order/type depends on template)
          */
         void timestep(ScalarType& rVar, ScalarType& nTerms);
         
      protected:

         /**
          * @brief Compute an iteration of scheme (actuall computation depends on type)
          */
         void doIteration(ScalarType& rVar, ScalarType& nTerms);

      private:
   };

   template <typename TSim, template <typename> class TMethod> ETDScheme<TSim, TMethod>::ETDScheme(DynamoFloat a, DynamoFloat b, const typename ETDScheme<TSim, TMethod>::BasisType &basis, TimestepParameters &tsteps, SmartSTrunc pTrunc)
      : TMethod<TSim>(a, b, basis, tsteps, pTrunc) 
   {
   }

   template <typename TSim, template <typename> class TMethod> ETDScheme<TSim, TMethod>::~ETDScheme()
   {
   }

   template <typename TSim, template <typename> class TMethod> void ETDScheme<TSim, TMethod>::timestep(typename ETDScheme<TSim, TMethod>::ScalarType& rVar, typename ETDScheme<TSim, TMethod>::ScalarType& nTerms)
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
         // Reorder the previous values to include newest step
         this->reorderPrevious(rVar, nTerms);

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

   template <typename TSim, template <typename> class TMethod> void ETDScheme<TSim, TMethod>::doIteration(typename ETDScheme<TSim, TMethod>::ScalarType& rVar, typename ETDScheme<TSim, TMethod>::ScalarType& nTerms)
   {
      if(this->hasIntermediate())
      {
         this->computeIntermediate(rVar, nTerms);
      } else
      {
         this->doStep(rVar, nTerms);
      }
   }

}

#endif // ETDSCHEME_HPP
