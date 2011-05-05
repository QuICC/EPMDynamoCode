/** \file ThetaRHSTOperatorSet.hpp
 *  \brief Implementation of the \f$\theta\f$ method RHS operator set
 */

#ifndef THETARHSTOPERATORSET_HPP
#define THETARHSTOPERATORSET_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "Domain/Truncation.hpp"
#include "Config/SimulationConfig.hpp"
#include "Timestepping/PredictorCorrector/Theta/ThetaTOperatorSet.hpp"

namespace EPMDynamo {

   /**
    * \brief Implementation of the \f$\theta\f$ method RHS operator set
    *
    * \tparam TSimType Type of the simulation
    */
   template <typename TSimType> class ThetaRHSTOperatorSet: public ThetaTOperatorSet<TSimType, SimulationConfig::OperatorType>
   {
      public:
         /// Typedef from Simulation trait to local radial basis type
         typedef typename TSimType::RadialBasisType    BasisType;

         /// Typedef from Simulation trait to local scalar type
         typedef typename TSimType::ScalarType    ScalarType;

         /**
          * @brief Constructor
          *
          * @param a Coefficient \f$a\f$ of timestep scheme
          * @param b Coefficient \f$b\f$ of timestep scheme
          * @param basis Reference to the basis used for the operators
          * @param pTrunc Truncation information
          * @param hasL0 Is the l=0 mode required?
          */
         ThetaRHSTOperatorSet(EPMFloat a, EPMFloat b, const BasisType &basis, SmartTruncation pTrunc, bool hasL0);

         /**
          * @brief Simple empty destructor
          */
         virtual ~ThetaRHSTOperatorSet() {};

         /**
          * @brief Update the timestep matrices
          *
          * @param dt   New timestep value
          */
         virtual void update(EPMFloat dt);
         
      protected:

      private:
   };

   template <typename TSimType> ThetaRHSTOperatorSet<TSimType>::ThetaRHSTOperatorSet(EPMFloat a, EPMFloat b, const typename ThetaRHSTOperatorSet<TSimType>::BasisType &basis, SmartTruncation pTrunc, bool hasL0)
      : ThetaTOperatorSet<TSimType, SimulationConfig::OperatorType>(a, b, basis, pTrunc, hasL0)
   {
   }

   template <typename TSimType> inline void ThetaRHSTOperatorSet<TSimType>::update(const EPMFloat dt)
   {
      // Set multiplicative factor for Laplacian
      EPMFloat factor = this->mB*(1.0-ThetaTraits<TSimType>::theta);

      // Set multiplicative factor for time matrix
      EPMFloat timeDiff = this->mA/dt;

      // Update stored operators
      this->updateOperators(factor, timeDiff);
   }

}

#endif // THETARHSTOPERATORSET_HPP
