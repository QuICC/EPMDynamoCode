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
          * @param pTrunc Truncation information
          */
         ETD3Operators(SmartTruncation pTrunc);

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
          *
          * \epmBug Wrong computation
          */
         virtual void update(const EPMFloat dt, const BasisType &basis);

         /**
          * @brief Update the ETD3 operators
          *
          * @param c Maximum eigen value
          * @param basis Radial basis 
          *
          * \epmBug Need to implement computation of inverse
          */
         void createOperators(const EPMFloat c, const BasisType &basis);

      protected:

      private:
   };

   template <typename TSimType> ETD3Operators<TSimType>::ETD3Operators(SmartTruncation pTrunc)
      : ETDNOperators<TSimType, 4>(pTrunc)
   {
   }

   template <typename TSimType> void ETD3Operators<TSimType>::update(const EPMFloat dt, const typename ETD0Operators<TSimType>::BasisType &basis)
   {
      this->createOperators(dt, basis);
   }

   template <typename TSimType> void ETD3Operators<TSimType>::createOperators(const EPMFloat c, const typename ETD3Operators<TSimType>::BasisType &basis)
   {
      // Storage for a temporary operator
      Matrix tmpM;
      Matrix tmpMInv;

      // Loop over all degrees
      for(int i = 0; i < this->etdF(0).nOp(); ++i)
      {
         // Define homogeneous operator
         this->rEtdF(0).rHarmOp(i).constructBOperator(c, basis.at(i).specLaplacian());

         // Store the operator including boundary conditions and its inverse
         tmpM = this->etdF(0).harmOp(i).op();
         //tmpMInv = INVERSE(tmpM);
         
         // Compute the exponential of the created operator
         this->computeScaledF0();

         // Initialise F1 numerator part of operator
         this->rEtdF(1).rHarmOp(i).rOp() = this->etdF(0).harmOp(i).op();

         // Remove identity
         this->rEtdF(1).rHarmOp(i).rOp().diagonal().cwise() -= 1.0;

         // initialise F2 numerator part of operator
         this->rEtdF(2).rHarmOp(i).rOp() = this->etdF(1).harmOp(i).op();

         // Remove c L
         this->rEtdF(2).rHarmOp(i).rOp() -= tmpM;

         // initialise F3 numerator part of operator
         this->rEtdF(3).rHarmOp(i).rOp() = this->etdF(2).harmOp(i).op();

         // Remove c^2 L^2/2
         tmpM *= tmpM;
         this->rEtdF(3).rHarmOp(i).rOp() -= tmpM.cwise()/2.0;

         // Multiply F1 by 1/c L^-1
         this->rEtdF(1).rHarmOp(i).rOp() *= tmpMInv;

         // Multiply F2 by 1/c^2 L^-2
         tmpMInv *= tmpMInv;
         this->rEtdF(2).rHarmOp(i).rOp() *= tmpMInv;

         // Multiply F3 by 1/c^3 L^-3
         tmpMInv *= tmpMInv;
         this->rEtdF(3).rHarmOp(i).rOp() *= tmpMInv;

         // Compute the unscale values
         this->computeSquaredF3();

         // Do finalisation step (for example factorisation)
         this->rEtdF(0).rHarmOp(i).finaliseOp();
         this->rEtdF(1).rHarmOp(i).finaliseOp();
         this->rEtdF(2).rHarmOp(i).finaliseOp();
         this->rEtdF(3).rHarmOp(i).finaliseOp();
      }
   }

}

#endif // ETD3OPERATORS_HPP
