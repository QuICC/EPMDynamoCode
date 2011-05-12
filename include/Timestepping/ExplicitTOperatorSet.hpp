/** \file ExplicitTOperatorSet.hpp
 *  \brief Implementation of the general implicit diffusion operator timestep operator set
 */

#ifndef EXPLICITTOPERATORSET_HPP
#define EXPLICITTOPERATORSET_HPP

// Configuration includes
//
#include "Config/SimulationConfig.hpp"

// System includes
//

// External includes
//

// Project includes
//
#include "General/EPMTypedefs.hpp"
#include "Operators/BoundedOperatorSet.hpp"
#include "Operators/BCRowOperator.hpp"

namespace EPMDynamo {

   /**
    * \brief Implementation of the general implicit diffusion operator timestep operator set
    *
    * \tparam TOpType Type of the operator
    */
   template <typename TOpType> class ExplicitTOperatorSet: public BoundedOperatorSet<TOpType, BCRowOperator>
   {
      public:
         /// Typedef from Simulation trait to local truncation type
         typedef SimulationConfig::NumericalScheme::RadialBasisType  BasisType;

         /**
          * @brief Constructor
          *
          * @param pTrunc Truncation information
          * @param hasL0 Is the l=0 mode required?
          */
         ExplicitTOperatorSet(SmartTruncation pTrunc, bool hasL0);

         /**
          * @brief Simple empty destructor
          */
         virtual ~ExplicitTOperatorSet() {};
         
      protected:
         /**
          * @brief Update the Real operator
          *
          * @param factor Multiplicative factor
          * @param timeDiff Timestep length
          * @param basis Radial basis
          */
         void createOperators(const EPMFloat factor, const EPMFloat timeDiff, const BasisType &basis);

      private:
   };

   template <typename TOpType> ExplicitTOperatorSet<TOpType>::ExplicitTOperatorSet(SmartTruncation pTrunc, bool hasL0)
      : BoundedOperatorSet<TOpType, BCRowOperator>(pTrunc, hasL0)
   {
   }

   template <typename TOpType> void ExplicitTOperatorSet<TOpType>::createOperators(const EPMFloat factor, const EPMFloat timeDiff, const typename ExplicitTOperatorSet<TOpType>::BasisType &basis)
   {
      // Storage for "id" of operator
      int l;

      // Loop over all degrees
      for(int i = this->minL(); i < this->nOp(); ++i)
      {
         // Get degree of the current operator
         l = this->harmOp(i).id();

         // Set the Laplacian component
         this->rHarmOp(i).rOp() = factor * basis.at(l).specLaplacian().op();

         // Add time derivative to diagonal
         this->rHarmOp(i).rOp().diagonal().array() += timeDiff;

         // Do finalisation step (for example factorisation)
         this->rHarmOp(i).finaliseOp();
      }
   }
}

#endif // EXPLICITTOPERATORSET_HPP
