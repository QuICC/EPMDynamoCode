/** \file PhysicalScalarBase.cpp
 *  \brief Implementation of physical scalar base
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "PhysicalFields/Fields/PhysicalScalarBase.hpp"

// Project includes
//

namespace EPMDynamo {
   
   PhysicalScalarBase::PhysicalScalarBase(SmartTruncation pTrunc, PhysicalScalarBase::TransformType &transform)
      : PhysicalRTPScalar(pTrunc, transform), mPerturbation(pTrunc, true)
   {
   }

   void PhysicalScalarBase::initialiseZeros()
   {
      PhysicalRTPScalar::initialiseZeros();

      this->mPerturbation.initialiseZeros();
   }

}
