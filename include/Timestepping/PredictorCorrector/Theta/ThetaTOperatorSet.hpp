/** \file ThetaTOperatorSet.hpp
 *  \brief General implementation of the \f$\theta\f$-method operator set
 */

#ifndef THETATOPERATORSET_HPP
#define THETATOPERATORSET_HPP

// Configuration includes
//
#include "Config/SimulationConfig.hpp"

// System includes
//
#include <vector>

// External includes
//

// Project includes
//
#include "General/EPMTypedefs.hpp"

namespace EPMDynamo {

   /**
    * \brief General implementation of the \f$\theta\f$-method operator set
    */
   template <typename TOpType, template <typename> class TSetType> class ThetaTOperatorSet: public TSetType<TOpType>
   {
      public:
         /// Typedef from Simulation trait to local truncation type
         typedef SimulationConfig::NumericalScheme::RadialBasisType  BasisType;

         /**
          * @brief Constructor
          *
          * @param a Coefficient \f$a\f$ of timestep scheme
          * @param b Coefficient \f$b\f$ of timestep scheme
          * @param basis Reference to the basis used for the operators
          * @param pTrunc Truncation information
          * @param hasL0 Is l=0 mode required?
          */
         ThetaTOperatorSet(EPMFloat a, EPMFloat b, const BasisType &basis, SmartTruncation pTrunc, bool hasL0);

         /**
          * @brief Simple empty destructor
          */
         virtual ~ThetaTOperatorSet() {};

         /**
          * @brief Update the timestep matrices
          *
          * Still a pure virtual function
          *
          * @param dt   New timestep value
          */
         virtual void update(const EPMFloat dt) = 0;

         /**
          * @brief Init the operators
          */
         void initOperators();
         
      protected:
         /**
          * @brief \f$a\f$ coefficient of timestep operator
          */
         EPMFloat    mA;

         /**
          * @brief \f$b\f$ coefficient of timestep operator
          */
         EPMFloat    mB;

         /**
          * @brief Reference to the radial basis
          */
         const BasisType&   mrBasis;

         /**
          * @brief Update the Real operator
          *
          * @param factor Multiplicative factor
          * @param timeDiff Timestep length
          */
         void updateOperators(const EPMFloat factor, const EPMFloat timeDiff);

      private:
         /**
          * @brief Lock for blocking setting of the theta parameter
          */
         bool mIsThetaLocked;
   };

   template <typename TOpType, template <typename> class TSetType> ThetaTOperatorSet<TOpType, TSetType>::ThetaTOperatorSet(EPMFloat a, EPMFloat b, const typename ThetaTOperatorSet<TOpType, TSetType>::BasisType &basis, SmartTruncation pTrunc, bool hasL0)
      : TSetType<TOpType>(pTrunc, hasL0), mA(a), mB(b), mrBasis(basis), mIsThetaLocked(false)
   {
   }

   template<typename TOpType, template <typename> class TSetType> void ThetaTOperatorSet<TOpType, TSetType>::initOperators()
   {
      // Forbid any further call to setTheta!
      this->mIsThetaLocked = true;

      // Call inherited initOperators method
      TSetType<TOpType>::initOperators();
   }

   template <typename TOpType, template <typename> class TSetType> void ThetaTOperatorSet<TOpType, TSetType>::updateOperators(const EPMFloat factor, const EPMFloat timeDiff)
   {
      this->createOperators(factor, timeDiff, this->mrBasis);
   }
}

#endif // THETATOPERATORSET_HPP
