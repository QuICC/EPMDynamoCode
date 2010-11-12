/** \file ImplicitTOperatorSet.hpp
 *  \brief Implementation of the general implicit diffusion operator timestep operator set
 */

#ifndef IMPLICITTOPERATORSET_HPP
#define IMPLICITTOPERATORSET_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "Simulations/Traits/SimulationTraits.hpp"
#include "General/EPMTypedefs.hpp"
#include "Operators/BoundedOperatorSet.hpp"

namespace EPMDynamo {

   /**
    * \brief Implementation of the general implicit diffusion operator timestep operator set
    *
    * \tparam TSimType Type of the simulation
    * \tparam TOpType Type of the operator
    */
   template <typename TSimType, typename TOpType> class ImplicitTOperatorSet: public BoundedOperatorSet<TOpType>
   {
      public:
         /// Typedef from Simulation trait to local truncation type
         typedef typename TSimType::RadialBasisType  BasisType;

         /**
          * @brief Constructor
          *
          * @param pTrunc Truncation information
          * @param hasL0 Is the l=0 mode required?
          */
         ImplicitTOperatorSet(SmartTruncation pTrunc, bool hasL0);

         /**
          * @brief Simple empty destructor
          */
         virtual ~ImplicitTOperatorSet() {};
         
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

   template <typename TSimType, typename TOpType> ImplicitTOperatorSet<TSimType, TOpType>::ImplicitTOperatorSet(SmartTruncation pTrunc, bool hasL0)
      : BoundedOperatorSet<TOpType>(pTrunc, hasL0)
   {
   }

   template <typename TSimType, typename TOpType> void ImplicitTOperatorSet<TSimType, TOpType>::createOperators(const EPMFloat factor, const EPMFloat timeDiff, const typename ImplicitTOperatorSet<TSimType, TOpType>::BasisType &basis)
   {
      // Storage for "id" of operator
      int l;

      // Loop over all degrees
      for(int i = this->minL(); i < this->nOp(); ++i)
      {
         // Get degree of the current operator
         l = this->harmOp(i).id();

         // Construct the bounded laplacian operator
         this->rHarmOp(i).constructBOperator(factor, basis.at(l).specLaplacian());

         // Add time derivative to diagonal
         this->rHarmOp(i).rOp().diagonal().array() += timeDiff;

         // Do finalisation step (for example factorisation)
         this->rHarmOp(i).finaliseOp();
      }
   }
}

#endif // IMPLICITTOPERATORSET_HPP
