/** \file IterativeScheme.hpp
 *  \brief Implementation of the high level representation of an iterative timestep scheme
 *
 *  This should allow a general implementation for RK, PC and ETD methods
 */

#ifndef ITERATIVESCHEME_HPP
#define ITERATIVESCHEME_HPP

// Configuration includes
//
#include "Config/SimulationConfig.hpp"

// System includes
//

// External includes
//

// Project includes
//
#include "Domain/Truncation.hpp"
#include "Timestepping/TimestepParameters.hpp"

namespace EPMDynamo {

   /**
    * \brief Implementation of the high level representation of an iterative timestep scheme
    *
    * \tparam TMethod The actual timestep method in use
    */
   template <typename TMethod> class IterativeScheme: public TMethod
   {
      public:
         /// Typedef from Simulation trait to local radial basis type
         typedef SimulationConfig::NumericalScheme::RadialBasisType    BasisType;

         /// Typedef from Simulation trait to local scalar type
         typedef SimulationConfig::NumericalScheme::ScalarType    ScalarType;

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
         IterativeScheme(EPMFloat a, EPMFloat b, const BasisType &basis, TimestepParameters &tsteps, SmartTruncation pTrunc, bool hasL0);

         /**
          * @brief Destructor
          */
         virtual ~IterativeScheme() {};

         /**
          * @brief Perform a time step (the actual order/type depends on the selected method)
          *
          * @param rVar Input/Output variable to timestep
          * @param nTerms Non linear terms
          */
         void timestep(ScalarType& rVar, ScalarType& nTerms);
         
      protected:

      private:
   };

   template <typename TMethod> IterativeScheme<TMethod>::IterativeScheme(EPMFloat a, EPMFloat b, const typename IterativeScheme<TMethod>::BasisType &basis, TimestepParameters &tsteps, SmartTruncation pTrunc, bool hasL0)
      : TMethod(a, b, basis, tsteps, pTrunc, hasL0) 
   {
   }

   template <typename TMethod> void IterativeScheme<TMethod>::timestep(typename IterativeScheme<TMethod>::ScalarType& rVar, typename IterativeScheme<TMethod>::ScalarType& nTerms)
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

#endif // ITERATIVESCHEME_HPP
