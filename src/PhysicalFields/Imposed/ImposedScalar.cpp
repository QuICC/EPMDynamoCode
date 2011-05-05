/** \file ImposedScalar.cpp
 *  \brief Implementation of scalar imposed field
 */

#ifndef IMPOSEDSCALAR_HPP
#define IMPOSEDSCALAR_HPP

// Configuration includes
//

// System includes
//

// External includes
//

// Project includes
//
#include "PhysicalFields/Imposed/ImposedScalar.hpp"

// Project includes
//

namespace EPMDynamo {

   ImposedScalar::ImposedScalar(SmartTruncation pTrunc, ImposedScalar::TransformType &transform)
      : PhysicalScalarBase(pTrunc, transform), mTotalField(pTrunc, true), mImposedField(pTrunc, true)
   {
   }

   void ImposedScalar::initialiseZeros()
   {
      PhysicalScalarBase::initialiseZeros();

      this->mImposedField.initialiseZeros();
   }

}
