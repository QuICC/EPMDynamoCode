/** \file PhysicalScalarGradientBase.cpp
 *  \brief Implementation of physical scalar base with gradient
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "PhysicalFields/Fields/PhysicalScalarGradientBase.hpp"

// Project includes
//

namespace EPMDynamo {
   
   PhysicalScalarGradientBase::PhysicalScalarGradientBase(SmartTruncation pTrunc, PhysicalScalarGradientBase::TransformType &transform)
      : PhysicalRTPScalarGradient(pTrunc, transform), mNeedGradTransform(0), mPerturbation(pTrunc, true)
   {
   }

   void PhysicalScalarGradientBase::initialiseZeros()
   {
      PhysicalRTPScalarGradient::initialiseZeros();

      this->mPerturbation.initialiseZeros();
   }

}
