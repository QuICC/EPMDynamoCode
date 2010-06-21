/** \file LaplacianBOperatorSet.hpp
 *  \brief Implementation of a set of Laplacian operators
 */

#ifndef LAPLACIANBOPERATORSET_HPP
#define LAPLACIANBOPERATORSET_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "Simulations/Traits/SimulationTraits.hpp"
#include "General/EPMTypedefs.hpp"
#include "Domain/Truncation.hpp"
#include "Operators/BoundedOperatorSet.hpp"

namespace EPMDynamo {

   /**
    * \brief Implementation of a set of Laplacian operators
    *
    * \tparam TSimType Type of the simulation
    * \tparam TOpType Type of the operator
    */
   template <typename TSimType, typename TOpType> class LaplacianBOperatorSet: public BoundedOperatorSet<TOpType>
   {
      public:
         /// Typedef from Simulation trait to local truncation type
         typedef typename TSimType::RadialBasisType  BasisType;

         /**
          * @brief Constructor
          *
          * @param basis Reference to the basis used for the operators
          * @param pTrunc Truncation information
          */
         LaplacianBOperatorSet(const BasisType &basis, SmartTruncation pTrunc);

         /**
          * @brief Simple empty destructor
          */
         virtual ~LaplacianBOperatorSet() {};
         
      protected:
         /**
          * @brief Reference to the radial basis
          */
         const BasisType&   mrBasis;

         /**
          * @brief Update the operators
          *
          * @param factor Multiplicative factor
          */
         void createOperators(const EPMFloat factor);

      private:
   };

   template <typename TSimType, typename TOpType> LaplacianBOperatorSet<TSimType, TOpType>::LaplacianBOperatorSet(const typename LaplacianBOperatorSet<TSimType, TOpType>::BasisType &basis, SmartTruncation pTrunc)
      : BoundedOperatorSet<TOpType>(pTrunc), mrBasis(basis)
   {
   }

   template <typename TSimType, typename TOpType> void LaplacianBOperatorSet<TSimType, TOpType>::createOperators(const EPMFloat factor)
   {
      // Storage for "id" of operator
      int l;

      // Loop over all degrees
      for(int i = 0; i < this->nOp(); ++i)
      {
         // Get degree of the current operator
         l = this->harmOp(i).id();

         // Construct the bounded operator
         this->rHarmOp(i).constructBOperator(factor, this->mrBasis.at(l).specLaplacian());

         // Do finalisation step (for example factorisation)
         this->rHarmOp(i).finaliseOp();
      }
   }
}

#endif // LAPLACIANBOPERATORSET_HPP
