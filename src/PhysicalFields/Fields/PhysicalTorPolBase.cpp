/** \file PhysicalTorPolBase.cpp
 *  \brief Base of the implementation of Toroidal/Poloidal expanded field field
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Project includes
//
#include "PhysicalFields/Fields/PhysicalTorPolBase.hpp"

// Project includes
//

namespace EPMDynamo {

   PhysicalTorPolBase::PhysicalTorPolBase(SmartTruncation pTrunc, PhysicalTorPolBase::TransformType &transform)
      : PhysicalRTPField(pTrunc, transform), mPerturbation(pTrunc)
   {
   }

   void PhysicalTorPolBase::initialiseZeros()
   {
      PhysicalRTPField::initialiseZeros();

      this->mPerturbation.initialiseZeros();
   }

}
