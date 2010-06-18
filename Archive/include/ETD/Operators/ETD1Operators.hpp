/** \file ETD1Operators.hpp
 */

#ifndef ETD1OPERATORS_HPP
#define ETD1OPERATORS_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "Domain/Truncation.hpp"
#include "Timestepping/TimestepParameters.hpp"
#include "Timestepping/ETD/Operators/ETDNOperators.hpp"
#include "Timestepping/ETD/ETDOperators.hpp"
#include "Simulations/Traits/SimulationTraits.hpp"

namespace EPMDynamo {

   template <typename TSim> class ETD1Operators: public ETDNOperators<TSim, 2>
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
         ETD1Operators(SmartSTrunc pTrunc);

         /**
          * @brief Destructor
          */
         virtual ~ETD1Operators();

         /**
          * @brief Update the ETD1 operators
          */
         void createOperators(const EPMFloat h, const BasisType &basis);

      protected:

      private:
   };

   template <typename TSim> ETD1Operators<TSim>::ETD0Operators(SmartSTrunc pTrunc)
      : ETDOperators<TSim,2>(pSTrunc)
   {
   }

   template <typename TSim> ETD1Operators<TSim>::~ETD0Operators()
   {
   }

   template <typename TSim> void ETD1Operators<TSim>::createOperators(const EPMFloat h, const typename ETD0Operators<TSim>::BasisType &basis)
   {
      // Storage for a temporary operator
      Matrix tmpM;

      // Storage for "id" of operator
      int l;

      // Loop over all degrees
      for(int i = 0; i < this->etdF(0).nOp(); ++i)
      {
         // Get degree of the current operator
         l = this->rEtdF(0).op(i).id();

         // Define homogeneous operator
         this->rEtdF(0).rOp(i).constructBOperator(c, basis.at(l).specLaplacian());

         // Store the operator including boundary conditions and its inverse
         tmpM = this->etdF(0).op(i).op();
         //tmpM = INVERSE(tmpM);

         // Compute the exponential of the created operator
         this->computeScaledF0();

         // Initialise F1 operator
         this->rEtdF(1).rOp(i).rOp() = this->etdF(0).op(i).op();

         // Remove identity
         this->rEtdF(1).rOp(i).rOp().diagonal().cwise() -= 1.0;

         // Multiply F1 by 1/c L^-1
         this->rEtdF(1).rOp(i).rOp() *= tmpM;

         // Compute the unscale values
         this->computeSquaredF1();

         // Do finalisation step (for example factorisation)
         this->rEtdF(0).rOp(i).finaliseOp();
         this->rEtdF(1).rOp(i).finaliseOp();
      }
   }

}

#endif // ETD1OPERATORS_HPP
