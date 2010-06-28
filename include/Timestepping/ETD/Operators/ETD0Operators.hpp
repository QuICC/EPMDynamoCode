/** \file ETD0Operators.hpp
 *  \brief Implementation of the ETD0 operators
 */

#ifndef ETD0OPERATORS_HPP
#define ETD0OPERATORS_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "Domain/Truncation.hpp"
#include "Timestepping/TimestepParameters.hpp"
#include "Timestepping/ETD/Operators/ETDNOperators.hpp"
#include "Simulations/Traits/SimulationTraits.hpp"

namespace EPMDynamo {

   /**
    * \brief Implementation of the ETD 0 operators
    *
    * \tparam TSimType Type of the simulation
    */
   template <typename TSimType> class ETD0Operators: public ETDNOperators<TSimType, 1>
   {
      public:
         /// Typedef from Simulation trait to local radial basis type
         typedef typename TSimType::RadialBasisType    BasisType;

         /// Typedef from Simulation trait to local scalar type
         typedef typename TSimType::ScalarType    ScalarType;

         /**
          * @brief Constructor
          *
          * @param pTrunc Truncation information
          */
         ETD0Operators(SmartTruncation pTrunc);

         /**
          * @brief Destructor
          */
         virtual ~ETD0Operators() {};

         /**
          * @brief Update the ETD0 operators
          *
          * @param c Maximum eigen value
          * @param basis Radial basis 
          *
          * \bug Documentation problem
          */
         void createOperators(const EPMFloat c, const BasisType &basis);

      protected:

      private:
   };

   template <typename TSimType> ETD0Operators<TSimType>::ETD0Operators(SmartTruncation pTrunc)
      : ETDOperators<TSimType, 1>(pTrunc)
   {
   }

   template <typename TSimType> void ETD0Operators<TSimType>::createOperators(const EPMFloat c, const typename ETD0Operators<TSimType>::BasisType &basis)
   {
      // Update the required scaling power
      this->updateScalings(c);

      // Storage for "id" of operator
      int l;

      // Loop over all degrees
      for(int i = 0; i < this->etdF(0).nOp(); ++i)
      {
         // Get degree of the current operator
         l = this->rEtdF(0).op(i).id();

         // Define homogeneous operator
         this->rEtdF(0).rOp(i).constructBOperator(c, basis.at(l).specLaplacian());

         // Compute the scaled exponential of the created operator
         this->computeScaledF0();

         // Compute the exponential from the possibly scaled exponential operator
         this->computeSquaredF0();

         // Do finalisation step (for example factorisation)
         this->rEtdF(0).rOp(i).finaliseOp();
      }
   }

}

#endif // ETD0OPERATORS_HPP
