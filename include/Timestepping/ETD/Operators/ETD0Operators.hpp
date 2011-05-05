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
#include "Config/SimulationConfig.hpp"

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
          * @param c Stiffness constant
          * @param pTrunc Truncation information
          * @param hasL0 Is l=0 mode required?
          */
         ETD0Operators(EPMFloat c, SmartTruncation pTrunc, bool hasL0);

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
          */
         virtual void update(const EPMFloat dt, const BasisType &basis);

         /**
          * @brief Update the ETD0 operators
          *
          * @param h Timestep length
          * @param basis Radial basis 
          */
         void createOperators(const EPMFloat h, const BasisType &basis);

      protected:

      private:
   };

   template <typename TSimType> ETD0Operators<TSimType>::ETD0Operators(EPMFloat c, SmartTruncation pTrunc, bool hasL0)
      : ETDNOperators<TSimType, 1>(c, pTrunc, hasL0)
   {
   }

   template <typename TSimType> void ETD0Operators<TSimType>::update(const EPMFloat dt, const typename ETD0Operators<TSimType>::BasisType &basis)
   {
      this->createOperators(dt, basis);
   }

   template <typename TSimType> void ETD0Operators<TSimType>::createOperators(const EPMFloat h, const typename ETD0Operators<TSimType>::BasisType &basis)
   {
      // Update the required scaling power
      this->updateScalings(h);

      // Loop over all degrees
      for(int l = this->etdF(0).minL(); l < this->etdF(0).nOp(); ++l)
      {
         // Define homogeneous operator
         this->rEtdF(0).rHarmOp(l).constructBOperator(h*this->c(), basis.at(l).specLaplacian());

         // Compute the scaled exponential of the created operator
         this->computeScaledF0();

         // Compute the exponential from the possibly scaled exponential operator
         this->computeSquaredF0();

         // Do finalisation step (for example factorisation)
         this->rEtdF(0).rHarmOp(l).finaliseOp();
      }
   }

}

#endif // ETD0OPERATORS_HPP
