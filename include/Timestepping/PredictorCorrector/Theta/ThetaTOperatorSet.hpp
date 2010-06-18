/** \file ThetaTOperatorSet.hpp
 *  \brief General implementation of the \f$\theta\f$-method operator set
 */

#ifndef THETATOPERATORSET_HPP
#define THETATOPERATORSET_HPP

// System includes
//
#include <vector>

// External includes
//

// Project includes
//
#include "Simulations/Traits/SimulationTraits.hpp"
#include "General/EPMTypedefs.hpp"
#include "General/EPMException.hpp"
#include "Timestepping/ImplicitTOperatorSet.hpp"

namespace EPMDynamo {

   /**
    * \brief General implementation of the \f$\theta\f$-method operator set
    */
   template <typename TSimType, typename TOpType> class ThetaTOperatorSet: public ImplicitTOperatorSet<TSimType, TOpType>
   {
      public:
         /// Typedef from Simulation trait to local truncation type
         typedef typename TSimType::RadialBasisType  BasisType;

         /**
          * @brief Constructor
          *
          * @param a Coefficient \f$a\f$ of timestep scheme
          * @param b Coefficient \f$b\f$ of timestep scheme
          * @param basis Reference to the basis used for the operators
          * @param pTrunc Truncation information
          */
         ThetaTOperatorSet(EPMFloat a, EPMFloat b, const BasisType &basis, SmartTruncation pTrunc);

         /**
          * @brief Simple empty destructor
          */
         virtual ~ThetaTOperatorSet() {};

         /**
          * @brief Get the theta parameter
          */
         EPMFloat theta() const;

         /**
          * @brief Update the timestep matrices
          *
          * Still a pure virtual function
          *
          * @param dt   New timestep value
          */
         virtual void update(const EPMFloat dt) = 0;

         /**
          * @brief Set the theta parameter of the scheme
          *
          * @param theta Implicitness
          *
          * \bug This has to be reimplemented
          */
         void setTheta(EPMFloat theta);

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
          * @brief Theta parameter of numerical scheme
          */
         EPMFloat   mTheta;

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

         /**
          * @brief Default theta parameter for all operators
          */
         static const EPMFloat   DEFAULT_THETA;

      private:
         /**
          * @brief Lock for blocking setting of the theta parameter
          */
         bool mIsThetaLocked;
   };

   template <typename TSimType, typename TOpType> ThetaTOperatorSet<TSimType, TOpType>::ThetaTOperatorSet(EPMFloat a, EPMFloat b, const typename ThetaTOperatorSet<TSimType, TOpType>::BasisType &basis, SmartTruncation pTrunc)
      : ImplicitTOperatorSet<TSimType, TOpType>(pTrunc), mA(a), mB(b), mTheta(DEFAULT_THETA), mrBasis(basis), mIsThetaLocked(false)
   {
   }
   
   template <typename TSimType, typename TOpType> inline EPMFloat ThetaTOperatorSet<TSimType, TOpType>::theta() const
   {
      return this->mTheta;
   }

   template<typename TSimType, typename TOpType> const EPMFloat ThetaTOperatorSet<TSimType, TOpType>::DEFAULT_THETA = 0.5;

   template<typename TSimType, typename TOpType> void ThetaTOperatorSet<TSimType, TOpType>::setTheta(EPMFloat theta)
   {
      // Check if simulation already started, the theta parameter should not be changed during run
      if(this->mIsThetaLocked)
      {
         throw EPMException("ThetaTOperatorSet::setTheta", "Tried to change theta after initialisation");
      } else
      {
         this->mTheta = theta;
      }
   }

   template<typename TSimType, typename TOpType> void ThetaTOperatorSet<TSimType, TOpType>::initOperators()
   {
      // Forbid any further call to setTheta!
      this->mIsThetaLocked = true;

      // Call inherited initOperators method
      ImplicitTOperatorSet<TSimType, TOpType>::initOperators();
   }

   template <typename TSimType, typename TOpType> void ThetaTOperatorSet<TSimType, TOpType>::updateOperators(const EPMFloat factor, const EPMFloat timeDiff)
   {
      this->createOperators(factor, timeDiff, this->mrBasis);
   }
}

#endif // THETATOPERATORSET_HPP
