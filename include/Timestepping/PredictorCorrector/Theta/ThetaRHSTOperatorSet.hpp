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
#include "Simulations/SimulationTraits.hpp"
#include "Timestepping/PredictorCorrector/Theta/ThetaTOperatorSet.hpp"

namespace EPMDynamo {

   /**
    * \brief Implementation of the \f$\theta\f$ method RHS operator set
    *
    * \tparam TSimType Type of the simulation
    *
    * \bug Needs to be reviewed and cleaned
    */
   template <typename TSimType> class ThetaRHSTOperatorSet: public ThetaTOperatorSet<TSimType, typename SimulationTraits<TSimType>::OperatorType>
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
          */
         ThetaRHSTOperatorSet(DynamoFloat a, DynamoFloat b, const BasisType &basis, SmartTruncation pTrunc);

         /**
          * @brief Simple empty destructor
          */
         virtual ~ThetaRHSTOperatorSet() {};

         /**
          * @brief Update the timestep matrices
          *
          * @param dt   New timestep value
          */
         virtual void update(DynamoFloat dt);
         
      protected:

      private:
   };

   template <typename TSimType> ThetaRHSTOperatorSet<TSimType>::ThetaRHSTOperatorSet(DynamoFloat a, DynamoFloat b, const typename ThetaRHSTOperatorSet<TSimType>::BasisType &basis, SmartTruncation pTrunc)
      : ThetaTOperatorSet<TSimType, typename SimulationTraits<TSimType>::OperatorType>(a, b, basis, pTrunc)
   {
   }

   template <typename TSimType> inline void ThetaRHSTOperatorSet<TSimType>::update(const DynamoFloat dt)
   {
      // Set multiplicative factor for Laplacian
      DynamoFloat factor = this->mB*(1.0-this->mTheta);

      // Set multiplicative factor for time matrix
      DynamoFloat timeDiff = this->mA/dt;

      // Update stored operators
      this->updateOperators(factor, timeDiff);
   }

}

#endif // THETARHSTOPERATORSET_HPP
