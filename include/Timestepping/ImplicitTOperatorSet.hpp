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
#include "Simulations/SimulationTraits.hpp"
#include "General/EPMTypedefs.hpp"
#include "Operators/BoundedOperatorSet.hpp"

namespace EPMDynamo {

   /**
    * \brief Implementation of the general implicit diffusion operator timestep operator set
    *
    * \tparam TSimType Type of the simulation
    * \tparam TOpType Type of the operator
    * \tparam TBCType Type of the boundary condition
    *
    * \bug Review this implementation regarding new operator implementation
    */
   template <typename TSimType, typename TOpType, BCType TBCType> class ImplicitTOperatorSet: public BoundedOperatorSet<TOpType, TBCType>
   {
      public:
         /// Typedef from Simulation trait to local truncation type
         typedef typename TSimType::RadialBasisType  BasisType;

         /**
          * @brief Constructor
          *
          */
         ImplicitTOperatorSet(SmartTruncation pTrunc);

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
         void createOperators(const DynamoFloat factor, const DynamoFloat timeDiff, const BasisType &basis);

      private:
   };

   template <typename TSimType, typename TOpType, BCType TBCType> ImplicitTOperatorSet<TSimType, TOpType, TBCType>::ImplicitTOperatorSet(SmartTruncation pTrunc)
      : BoundedOperatorSet<TOpType, TBCType>(pTrunc)
   {
   }

   template <typename TSimType, typename TOpType, BCType TBCType> void ImplicitTOperatorSet<TSimType, TOpType, TBCType>::createOperators(const DynamoFloat factor, const DynamoFloat timeDiff, const typename ImplicitTOperatorSet<TSimType, TOpType, TBCType>::BasisType &basis)
   {
      // Storage for "id" of operator
      int l;

      // Loop over all degrees
      for(int i = 0; i < this->nOp(); ++i)
      {
         // Get degree of the current operator
         l = this->op(i).id();

         // Construct the bounded laplacian operator
         this->rOp(i).constructBOperator(factor, basis.at(l).specLaplacian());

         // Add time derivative to diagonal
         this->rOp(i).rOp().diagonal().cwise() += timeDiff;

         // Do finalisation step (for example factorisation)
         this->rOp(i).finaliseOp();
      }
   }
}

#endif // IMPLICITTOPERATORSET_HPP
