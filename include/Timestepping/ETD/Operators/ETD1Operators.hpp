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

   template <typename TSimType> ETD1Operators<TSimType>::ETD1Operators(EPMFloat c, SmartTruncation pTrunc, bool hasL0)
      : ETDNOperators<TSimType,2>(c, pTrunc, hasL0)
   {
   }

   template <typename TSimType> void ETD1Operators<TSimType>::update(const EPMFloat dt, const typename ETD1Operators<TSimType>::BasisType &basis)
   {
      this->createOperators(dt, basis);
   }

   template <typename TSimType> void ETD1Operators<TSimType>::createOperators(const EPMFloat h, const typename ETD1Operators<TSimType>::BasisType &basis)
   {
      // Update the required scaling power
      this->updateScalings(h);

      // Loop over all degrees
      for(int l = this->etdF(0).minL(); l < this->etdF(0).nOp(); ++l)
      {
         // Define homogeneous operator
         this->rEtdF(0).rHarmOp(l).constructBOperator(h*this->c(), basis.at(l).specLaplacian());

         // Compute the taylor expansion of the Fk functionals of the created operator
         this->computeScaledF1(l);

         // Compute the unscaled values
         this->computeSquaredF1(l);

         // Include the missing h factor
         this->rEtdF(1).rHarmOp(l).rOp() *= h;

         // Correct the results for zero eigenvalues
         if(!this->isFullRank(l))
         { 
            this->rEtdF(1).rHarmOp(l).rOp().topRows(1) += h*this->etdF(0).harmOp(l).op().leftCols(1).transpose();
         }

         // Do finalisation step
         this->rEtdF(0).rHarmOp(l).finaliseOp();
         this->rEtdF(1).rHarmOp(l).finaliseOp();
      }
   }

}

#endif // ETD1OPERATORS_HPP
