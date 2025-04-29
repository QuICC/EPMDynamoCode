/** \file PhysicalBase.cpp
 *  \brief Base of the implementation of the physical fields
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "PhysicalFields/Fields/PhysicalBase.hpp"

// Project includes
//

namespace EPMDynamo {

   PhysicalBase::PhysicalBase(SmartTruncation pTrunc, PhysicalBase::TransformType &transform)
      : mrTransform(transform), mNeedTransform(0), mpTrunc(pTrunc)
   {
   }

}
