/** \file ETD3Operators.hpp
 *  \brief Implementation of the ETD3 operators
 */

#ifndef ETD3OPERATORS_HPP
#define ETD3OPERATORS_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "Domain/Truncation.hpp"
#include "Timestepping/TimestepParameters.hpp"
#include "Timestepping/ETD/ETDOperators.hpp"
#include "Timestepping/ETD/ETDSchemeTraits.hpp"
#include "Timestepping/ETD/Operators/ETDNOperators.hpp"
#include "Simulations/Traits/SimulationTraits.hpp"

namespace EPMDynamo {

   /**
    * \brief Implementation of the ETD3 operators
    *
    * \tparam TSimType Type of the simulation
    */
   template <typename TSimType> class ETD3Operators: public ETDNOperators<TSimType, 4>
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
         ETD3Operators(EPMFloat c, SmartTruncation pTrunc, bool hasL0);

         /**
          * @brief Destructor
          */
         virtual ~ETD3Operators() {};

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
          * @brief Update the ETD3 operators
          *
          * @param h Timestep length
          * @param basis Radial basis 
          */
         void createOperators(const EPMFloat h, const BasisType &basis);

      protected:

      private:
   };

   template <typename TSimType> ETD3Operators<TSimType>::ETD3Operators(EPMFloat c, SmartTruncation pTrunc, bool hasL0)
      : ETDNOperators<TSimType, 4>(c, pTrunc, hasL0)
   {
   }

   template <typename TSimType> void ETD3Operators<TSimType>::update(const EPMFloat dt, const typename ETD0Operators<TSimType>::BasisType &basis)
   {
      this->createOperators(dt, basis);
   }

   template <typename TSimType> void ETD3Operators<TSimType>::createOperators(const EPMFloat h, const typename ETD3Operators<TSimType>::BasisType &basis)
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

         // Store the operator including boundary conditions and its inverse
         tmpM = this->etdF(0).harmOp(i).op();

         // Scale the Lh operator
         this->scaleOperator(tmpM, l);
         
         // Compute the exponential of the created operator
         this->computeScaledF0();

         // Initialise F1 numerator part of operator
         this->rEtdF(1).rHarmOp(i).rOp() = this->etdF(0).harmOp(i).op();

         // Remove identity
         this->rEtdF(1).rHarmOp(i).rOp().diagonal().array() -= 1.0;

         // initialise F2 numerator part of operator
         this->rEtdF(2).rHarmOp(i).rOp() = this->etdF(1).harmOp(i).op();

         // Remove c L
         this->rEtdF(2).rHarmOp(i).rOp() -= tmpM;

         // initialise F3 numerator part of operator
         this->rEtdF(3).rHarmOp(i).rOp() = this->etdF(2).harmOp(i).op();

         // Remove c^2 L^2/2
         this->rEtdF(3).rHarmOp(i).rOp() -= (tmpM*tmpM).array()/2.0;

         // Compute inverse of M
         tmpM *= h;
         this->computeInverse(tmpM, l);

         // Multiply F1 by 1/c L^-1
         this->rEtdF(1).rHarmOp(i).rOp() *= tmpM;

         // Multiply F2 by 1/c^2 L^-2
         tmpM *= tmpM;
         this->rEtdF(2).rHarmOp(i).rOp() *= tmpM;

         // Multiply F3 by 1/c^3 L^-3
         tmpM *= tmpM;
         this->rEtdF(3).rHarmOp(i).rOp() *= tmpM;

         // Compute the unscale values
         this->computeSquaredF3();

         // Include the missing h factor
         this->rEtdF(1).rHarmOp(l).rOp() *= h;

         // Include the missing h^2 factor
         this->rEtdF(2).rHarmOp(l).rOp() *= h*h;

         // Include the missing h^3 factor
         this->rEtdF(3).rHarmOp(l).rOp() *= h*h*h;

         // Do finalisation step (for example factorisation)
         this->rEtdF(0).rHarmOp(i).finaliseOp();
         this->rEtdF(1).rHarmOp(i).finaliseOp();
         this->rEtdF(2).rHarmOp(i).finaliseOp();
         this->rEtdF(3).rHarmOp(i).finaliseOp();
      }
   }

}

#endif // ETD3OPERATORS_HPP
