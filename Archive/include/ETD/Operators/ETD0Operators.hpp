/** \file ETD0Operators.hpp
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

   template <typename TSim> class ETD0Operators: public ETDNOperators<TSim, 1>
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
         ETD0Operators(SmartSTrunc pTrunc);

         /**
          * @brief Destructor
          */
         virtual ~ETD0Operators();

         /**
          * @brief Update the ETD0 operators
          */
         void createOperators(const EPMFloat c, const BasisType &basis);

      protected:

      private:
   };

   template <typename TSim> ETD0Operators<TSim>::ETD0Operators(SmartSTrunc pTrunc)
      : ETDOperators<TSim, 1>(pSTrunc)
   {
   }

   template <typename TSim> ETD0Operators<TSim>::~ETD0Operators()
   {
   }

   template <typename TSim> void ETD0Operators<TSim>::createOperators(const EPMFloat c, const typename ETD0Operators<TSim>::BasisType &basis)
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
