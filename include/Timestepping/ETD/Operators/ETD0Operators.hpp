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
          * @brief Update the timestep matrices
          *
          * Still a pure virtual function
          *
          * @param dt   New timestep value
          * @param basis Radial basis
          *
          * \epmBug Wrong computation, there is no conversion between dt and c
          */
         virtual void update(const EPMFloat dt, const BasisType &basis);

         /**
          * @brief Update the ETD0 operators
          *
          * @param c Maximum eigen value
          * @param basis Radial basis 
          */
         void createOperators(const EPMFloat c, const BasisType &basis);

      protected:

      private:
   };

   template <typename TSimType> ETD0Operators<TSimType>::ETD0Operators(SmartTruncation pTrunc)
      : ETDNOperators<TSimType, 1>(pTrunc)
   {
   }

   template <typename TSimType> void ETD0Operators<TSimType>::update(const EPMFloat dt, const typename ETD0Operators<TSimType>::BasisType &basis)
   {
      this->createOperators(dt, basis);
   }

   template <typename TSimType> void ETD0Operators<TSimType>::createOperators(const EPMFloat c, const typename ETD0Operators<TSimType>::BasisType &basis)
   {
      // Update the required scaling power
      this->updateScalings(c);

      // Loop over all degrees
      for(int i = 0; i < this->etdF(0).nOp(); ++i)
      {
         // Define homogeneous operator
         this->rEtdF(0).rHarmOp(i).constructBOperator(c, basis.at(i).specLaplacian());

         // Compute the scaled exponential of the created operator
         this->computeScaledF0();

         // Compute the exponential from the possibly scaled exponential operator
         this->computeSquaredF0();

         // Do finalisation step (for example factorisation)
         this->rEtdF(0).rHarmOp(i).finaliseOp();
      }
   }

}

#endif // ETD0OPERATORS_HPP
