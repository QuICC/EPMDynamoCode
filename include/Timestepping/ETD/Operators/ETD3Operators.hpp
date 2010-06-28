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

         /// Typedef for a smart pointer to ETDOperators
         typedef typename ETDOperators<TSimType> SmartETDOperators;

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
          * @brief Update the ETD3 operators
          *
          * @param h ??? Maximum eigen value ???
          * @param basis Radial basis 
          *
          * \epmBug Documentation problem
          */
         void createOperators(const EPMFloat h, const BasisType &basis);

      protected:

      private:
   };

   template <typename TSimType> ETD3Operators<TSimType>::ETD3Operators(SmartTruncation pTrunc)
      : ETDNOperators<TSimType, 4>(pTrunc)
   {
   }

   template <typename TSimType> void ETD3Operators<TSimType>::createOperators(const EPMFloat h, const typename ETD3Operators<TSimType>::BasisType &basis)
   {
      // Get list of degrees
      ArrayI   degrees = this->pIdxCtrl()->opLs();

      // Storage for a temporary operator
      Matrix tmpM;
      Matrix tmpMInv;

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
         //tmpMInv = INVERSE(tmpM);
         
         // Compute the exponential of the created operator
         this->computeScaledF0();

         // Initialise F1 numerator part of operator
         this->rEtdF(1).rOp(i).rOp() = this->etdF(0).op(i).op();

         // Remove identity
         this->rEtdF(1).rOp(i).rOp().diagonal().cwise() -= 1.0;

         // initialise F2 numerator part of operator
         this->rEtdF(2).rOp(i).rOp() = this->etdF(1).op(i).op();

         // Remove c L
         this->rEtdF(2).rOp(i).rOp() -= tmpM;

         // initialise F3 numerator part of operator
         this->rEtdF(3).rOp(i).rOp() = this->etdF(2).op(i).op();

         // Remove c^2 L^2/2
         tmpM *= tmpM;
         this->rEtdF(3).rOp(i).rOp() -= tmpM.cwise()/2.0;

         // Multiply F1 by 1/c L^-1
         this->rEtdF(1).rOp(i).rOp() *= tmpMInv;

         // Multiply F2 by 1/c^2 L^-2
         tmpMInv *= tmpMInv;
         this->rEtdF(2).rOp(i).rOp() *= tmpMInv;

         // Multiply F3 by 1/c^3 L^-3
         tmpMInv *= tmpMInv;
         this->rEtdF(3).rOp(i).rOp() *= tmpMInv;

         // Compute the unscale values
         this->computeSquaredF3();

         // Do finalisation step (for example factorisation)
         this->rEtdF(0).rOp(i).finaliseOp();
         this->rEtdF(1).rOp(i).finaliseOp();
         this->rEtdF(2).rOp(i).finaliseOp();
         this->rEtdF(3).rOp(i).finaliseOp();
      }
   }

}

#endif // ETD3OPERATORS_HPP
