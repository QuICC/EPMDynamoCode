/** \file LaplacianBOperatorSet.hpp
 *  \brief Implementation of a set of Laplacian operators
 */

#ifndef LAPLACIANBOPERATORSET_HPP
#define LAPLACIANBOPERATORSET_HPP

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
#include "Domain/Truncation.hpp"
#include "Operators/BoundedOperatorSet.hpp"

namespace EPMDynamo {

   /**
    * \brief Implementation of a set of Laplacian operators
    *
    * \tparam TOpType Type of the operator
    */
   template <typename TOpType> class LaplacianBOperatorSet: public BoundedOperatorSet<TOpType>
   {
      public:
         /// Typedef from Simulation trait to local truncation type
         typedef SimulationConfig::NumericalScheme::RadialBasisType  BasisType;

         /**
          * @brief Constructor
          *
          * @param basis Reference to the basis used for the operators
          * @param pTrunc Truncation information
          * @param hasL0 Is the l=0 mode required?
          */
         LaplacianBOperatorSet(const BasisType &basis, SmartTruncation pTrunc, bool hasL0);

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

   template <typename TOpType> LaplacianBOperatorSet<TOpType>::LaplacianBOperatorSet(const typename LaplacianBOperatorSet<TOpType>::BasisType &basis, SmartTruncation pTrunc, bool hasL0)
      : BoundedOperatorSet<TOpType>(pTrunc, hasL0), mrBasis(basis)
   {
   }

   template <typename TOpType> void LaplacianBOperatorSet<TOpType>::createOperators(const EPMFloat factor)
   {
      // Storage for "id" of operator
      int l;

      // Loop over all degrees
      for(int i = this->minL(); i < this->nOp(); ++i)
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
