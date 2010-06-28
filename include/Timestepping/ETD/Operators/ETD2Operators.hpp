/** \file ETD2Operators.hpp
 *  \brief Implementation of the ETD2 operators
 */

#ifndef ETD2OPERATORS_HPP
#define ETD2OPERATORS_HPP

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
    * \brief Implementation of the ETD2 operators
    *
    */
   template <typename TSimType> class ETD2Operators: public ETDNOperators<TSimType, 3>
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
         ETD2Operators(SmartTruncation pTrunc);

         /**
          * @brief Destructor
          */
         virtual ~ETD2Operators() {};

         /**
          * @brief Update the ETD2 operators
          *
          * @param c ?? Maximum eigen value ??
          * @param basis Radial basis 
          *
          * \epmBug Documentation problem
          */
         void createOperators(const EPMFloat c, const BasisType &basis);

      protected:

      private:
   };

   template <typename TSimType> ETD2Operators<TSimType>::ETD2Operators(SmartTruncation pTrunc)
      : ETDNOperators<TSimType, 3>(pTrunc)
   {
   }

   template <typename TSimType> void ETD2Operators<TSimType>::createOperators(const EPMFloat c, const typename ETD2Operators<TSimType>::BasisType &basis)
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

         // Compute inverse of M
         //tmpM = INVERSE(tmpM);

         // Multiply F1 by 1/c L^-1
         this->rEtdF(1).rOp(i).rOp() *= tmpM;

         // Multiply F2 by 1/c² L^-2
         tmpM *= tmpM;
         this->rEtdF(2).rOp(i).rOp() *= tmpM;

         // Compute the unscale values
         this->computeSquaredF2();

         // Do finalisation step (for example factorisation)
         this->rEtdF(0).rOp(i).finaliseOp();
         this->rEtdF(1).rOp(i).finaliseOp();
         this->rEtdF(2).rOp(i).finaliseOp();
      }
   }

}

#endif // ETD2OPERATORS_HPP
