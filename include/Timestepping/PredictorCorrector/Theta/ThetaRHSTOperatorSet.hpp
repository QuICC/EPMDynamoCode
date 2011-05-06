/** \file ThetaRHSTOperatorSet.hpp
 *  \brief Implementation of the \f$\theta\f$ method RHS operator set
 */

#ifndef THETARHSTOPERATORSET_HPP
#define THETARHSTOPERATORSET_HPP

// Configuration includes
//
#include "Config/SimulationConfig.hpp"
#include "Config/NumericalSchemeInc.hpp"

// System includes
//

// External includes
//

// Project includes
//
#include "Domain/Truncation.hpp"
#include "Timestepping/PredictorCorrector/Theta/ThetaTraits.hpp"
#include "Timestepping/PredictorCorrector/Theta/ThetaTOperatorSet.hpp"

namespace EPMDynamo {

   /**
    * \brief Implementation of the \f$\theta\f$ method RHS operator set
    */
   class ThetaRHSTOperatorSet: public ThetaTOperatorSet<SimulationConfig::OperatorType>
   {
      public:
         /// Typedef from Simulation trait to local radial basis type
         typedef SimulationConfig::NumericalScheme::RadialBasisType    BasisType;

         /// Typedef from Simulation trait to local scalar type
         typedef SimulationConfig::NumericalScheme::ScalarType    ScalarType;

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

   inline void ThetaRHSTOperatorSet::update(const EPMFloat dt)
   {
      // Set multiplicative factor for Laplacian
      EPMFloat factor = this->mB*(1.0-ThetaTraits::theta);

      // Set multiplicative factor for time matrix
      EPMFloat timeDiff = this->mA/dt;

      // Update stored operators
      this->updateOperators(factor, timeDiff);
   }

}

#endif // THETARHSTOPERATORSET_HPP
