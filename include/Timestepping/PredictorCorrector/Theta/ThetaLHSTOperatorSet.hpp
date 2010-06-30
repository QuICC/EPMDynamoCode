/** \file ThetaLHSTOperatorSet.hpp
 *  \brief Implementation of the \f$\theta\f$-method LHS operator set
 */

#ifndef THETALHSTOPERATORSET_HPP
#define THETALHSTOPERATORSET_HPP

// System includes
//
#include <vector>

// External includes
//

// Project includes
//
#include "Simulations/Traits/SimulationTraits.hpp"
#include "Timestepping/PredictorCorrector/Theta/ThetaTOperatorSet.hpp"
#include "Timestepping/PredictorCorrector/Theta/ThetaTraits.hpp"

namespace EPMDynamo {

   /**
    * \brief Implementation of the \f$\theta\f$-method LHS operator set
    *
    * \tparam TSimType Type of the simulation
    */
   template <typename TSimType> class ThetaLHSTOperatorSet: public ThetaTOperatorSet<TSimType, typename SimulationTraits<TSimType>::FactoredOpType>
   {
      public:
         /// Typedef from Simulation trait to local truncation type
         typedef typename TSimType::RadialBasisType    BasisType;

         /// Typedef from Simulation trait to local truncation type
         typedef typename TSimType::ScalarType    ScalarType;

         /**
          * @brief Constructor
          *
          * @param a Coefficient \f$a\f$ of timestep scheme
          * @param b Coefficient \f$b\f$ of timestep scheme
          * @param basis Reference to the basis used for the operators
          * @param pTrunc Truncation information
          * @param hasL0 Is the l=0 mode required ?
          */
         ThetaLHSTOperatorSet(EPMFloat a, EPMFloat b, const BasisType &basis, SmartTruncation pTrunc, bool hasL0);

         /**
          * @brief Simple empty destructor
          */
         virtual ~ThetaLHSTOperatorSet() {};

         /**
          * @brief Update the timestep matrices
          *
          * @param dt   New timestep value
          */
         virtual void update(EPMFloat dt);

         /**
          * @brief Solve timestep equation for unknown
          *
          * @param rVar Input/Output variable
          */
         void solve(ScalarType &rVar);

         /**
          * @brief Solve timestep equation for unkown using homogenous BC whatever real BC would be
          *
          * @param rVar Input/Output variable
          */
         void solveZero(ScalarType &rVar);
         
      protected:

      private:
   };

   template <typename TSimType> ThetaLHSTOperatorSet<TSimType>::ThetaLHSTOperatorSet(EPMFloat a, EPMFloat b, const typename ThetaLHSTOperatorSet<TSimType>::BasisType &basis, SmartTruncation pTrunc, bool hasL0)
      : ThetaTOperatorSet<TSimType, typename SimulationTraits<TSimType>::FactoredOpType>(a, b, basis, pTrunc, hasL0)
   {
   }

   template <typename TSimType> inline void ThetaLHSTOperatorSet<TSimType>::update(const EPMFloat dt)
   {
      // Set Laplacian multiplicative factor
      EPMFloat factor = -this->mB*ThetaTraits<TSimType>::theta;

      // Set time matrix multiplicative factor
      EPMFloat timeDiff = this->mA/dt;

      // Update stored operators
      this->updateOperators(factor, timeDiff);
   }

   template <typename TSimType> inline void ThetaLHSTOperatorSet<TSimType>::solve(typename ThetaLHSTOperatorSet<TSimType>::ScalarType &rVar)
   {
      // Get number of harmonic degrees
      int nL = rVar.nL();
      const int l0 = rVar.minL();

      // Loop over degrees to solve for all modes
      for(int l = l0; l <  nL; ++l)
      {
         this->solveOrders(rVar.rLShell(l), l);
      }
   }

   template <typename TSimType> inline void ThetaLHSTOperatorSet<TSimType>::solveZero(typename ThetaLHSTOperatorSet<TSimType>::ScalarType &rVar)
   {
      // Get number of harmonic degrees
      int nL = rVar.nL();
      const int l0 = rVar.minL();

      // Loop over degrees to solve for all modes with homogeneous operators
      for(int l = l0; l <  nL; ++l)
      {
         this->solveZeroOrders(rVar.rLShell(l), l);
      }
   }

}

#endif // THETALHSTOPERATORSET_HPP
