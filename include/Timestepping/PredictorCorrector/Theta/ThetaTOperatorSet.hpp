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
#include "Simulations/SimulationTraits.hpp"
#include "General/EPMTypedefs.hpp"
#include "General/EPMException.hpp"
#include "Timestepping/ImplicitTOperatorSet.hpp"

namespace EPMDynamo {

   /**
    * \brief General implementation of the \f$\theta\f$-method operator set
    */
   template <typename TSim, typename TOpType> class ThetaTOperatorSet: public ImplicitTOperatorSet<TSim, TOpType>
   {
      public:
         /// Typedef from Simulation trait to local truncation type
         typedef typename TSim::RadialBasisType  BasisType;

         /**
          * @brief Constructor
          *
          * @param a Coefficient \f$a\f$ of timestep scheme
          * @param b Coefficient \f$b\f$ of timestep scheme
          * @param basis Reference to the basis used for the operators
          * @param pTrunc Truncation information
          */
         ThetaTOperatorSet(DynamoFloat a, DynamoFloat b, const BasisType &basis, SmartTruncation pTrunc);

         /**
          * @brief Simple empty destructor
          */
         virtual ~ThetaTOperatorSet() {};

         /**
          * @brief Get the theta parameter
          */
         DynamoFloat theta() const;

         /**
          * @brief Update the timestep matrices
          *
          * Still a pure virtual function
          *
          * @param dt   New timestep value
          */
         virtual void update(const DynamoFloat dt) = 0;

         /**
          * @brief Set the theta parameter of the scheme
          *
          * @param theta Implicitness
          *
          * \bug This has to be reimplemented
          */
         void setTheta(DynamoFloat theta);

         /**
          * @brief Init the operators
          */
         void initOperators();
         
      protected:
         /**
          * @brief \f$a\f$ coefficient of timestep operator
          */
         DynamoFloat    mA;

         /**
          * @brief \f$b\f$ coefficient of timestep operator
          */
         DynamoFloat    mB;

         /**
          * @brief Theta parameter of numerical scheme
          */
         DynamoFloat   mTheta;

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
         void updateOperators(const DynamoFloat factor, const DynamoFloat timeDiff);

         /**
          * @brief Default theta parameter for all operators
          */
         static const DynamoFloat   DEFAULT_THETA;

      private:
         /**
          * @brief Lock for blocking setting of the theta parameter
          */
         bool mIsThetaLocked;
   };

   template <typename TSim, typename TOpType> ThetaTOperatorSet<TSim, TOpType>::ThetaTOperatorSet(DynamoFloat a, DynamoFloat b, const typename ThetaTOperatorSet<TSim, TOpType>::BasisType &basis, SmartTruncation pTrunc)
      : ImplicitTOperatorSet<TSim, TOpType>(pTrunc), mA(a), mB(b), mTheta(DEFAULT_THETA), mrBasis(basis), mIsThetaLocked(false)
   {
   }
   
   template <typename TSim, typename TOpType> inline DynamoFloat ThetaTOperatorSet<TSim, TOpType>::theta() const
   {
      return this->mTheta;
   }

   template<typename TSim, typename TOpType> const DynamoFloat ThetaTOperatorSet<TSim, TOpType>::DEFAULT_THETA = 0.5;

   template<typename TSim, typename TOpType> void ThetaTOperatorSet<TSim, TOpType>::setTheta(DynamoFloat theta)
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

   template<typename TSim, typename TOpType> void ThetaTOperatorSet<TSim, TOpType>::initOperators()
   {
      // Forbid any further call to setTheta!
      this->mIsThetaLocked = true;

      // Call inherited initOperators method
      ImplicitTOperatorSet<TSim, TOpType>::initOperators();
   }

   template <typename TSim, typename TOpType> void ThetaTOperatorSet<TSim, TOpType>::updateOperators(const DynamoFloat factor, const DynamoFloat timeDiff)
   {
      this->createOperators(factor, timeDiff, this->mrBasis);
   }
}

#endif // THETATOPERATORSET_HPP
