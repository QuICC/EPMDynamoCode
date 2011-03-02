/** \file ETD1Operators.hpp
 *  \brief Implementation of the ETD0 operators
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

   /**
    * \brief Implementation of the ETD 1 operators
    *
    * \tparam TSimType Type of the simulation
    */
   template <typename TSimType> class ETD1Operators: public ETDNOperators<TSimType, 2>
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
         ETD1Operators(EPMFloat c, SmartTruncation pTrunc, bool hasL0);

         /**
          * @brief Destructor
          */
         virtual ~ETD1Operators() {};

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
          * @brief Update the ETD1 operators
          *
          * @param h Timestep length
          * @param basis Radial basis 
          */
         void createOperators(const EPMFloat h, const BasisType &basis);

      protected:

      private:
   };

   template <typename TSimType> ETD1Operators<TSimType>::ETD0Operators(EPMFloat c, SmartTruncation pTrunc, bool hasL0)
      : ETDNOperators<TSimType,2>(c, pTrunc, hasL0)
   {
   }

   template <typename TSimType> void ETD1Operators<TSimType>::update(const EPMFloat dt, const typename ETD0Operators<TSimType>::BasisType &basis)
   {
      this->createOperators(dt, basis);
   }

   template <typename TSimType> void ETD1Operators<TSimType>::createOperators(const EPMFloat h, const typename ETD0Operators<TSimType>::BasisType &basis)
   {
      // Update the required scaling power
      this->updateScalings(h);

      // Storage for a temporary operator
      Matrix tmpM;

      // Loop over all degrees
      for(int i = this->etdF(0).minL(); i < this->etdF(0).nOp(); ++i)
      {
         // Define homogeneous operator
         this->rEtdF(0).rHarmOp(i).constructBOperator(h*this->c(), basis.at(i).specLaplacian());

         // Store the operator including boundary conditions
         tmpM = this->etdF(0).harmOp(i).op();

         // Scale the Lh operator
         this->scaleOperator(tmpM, l);

         // Compute the exponential of the created operator
         this->computeScaledF0();

         // Initialise F1 operator
         this->rEtdF(1).rHarmOp(i).rOp() = this->etdF(0).op(i).op();

         // Remove identity
         this->rEtdF(1).rHarmOp(i).rOp().diagonal().array() -= 1.0;

         // Compute inverse of operator
         tmpM *= h;
         this->computeInverse(tmpM);

         // Multiply F1 by 1/c L^-1
         this->rEtdF(1).rHarmOp(i).rOp() *= tmpM;

         // Compute the unscale values
         this->computeSquaredF1();

         // Include the missing h factor
         this->rEtdF(1).rHarmOp(l).rOp() *= h;

         // Do finalisation step (for example factorisation)
         this->rEtdF(0).rHarmOp(i).finaliseOp();
         this->rEtdF(1).rHarmOp(i).finaliseOp();
      }
   }

}

#endif // ETD1OPERATORS_HPP
