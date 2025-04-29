/** \file PhysicalQSTBase.cpp
 *  \brief Implementation of the base of the QST expanded physical field
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "PhysicalFields/Fields/PhysicalQSTBase.hpp"

// Project includes
//

namespace EPMDynamo {

   PhysicalQSTBase::PhysicalQSTBase(SmartTruncation pTrunc, PhysicalQSTBase::TransformType &transform)
      : PhysicalRTPField(pTrunc, transform), mPerturbation(pTrunc)
   {
   }

   void PhysicalQSTBase::initialiseZeros()
   {
      PhysicalRTPField::initialiseZeros();

      this->mPerturbation.initialiseZeros();
   }
}
