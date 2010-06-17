/** \file WorkflowBase.cpp
 *  \brief Source of base of implementation of workflow control
 */

// System includes
//

// External includes
//

// Class include
//
#include "Parallelisers/WorkflowBase.hpp"

// Project includes
//
#include "General/EPMException.hpp"

namespace EPMDynamo {

   WorkflowBase::WorkflowBase(const int nCore)
      : mNcore(nCore), mCoreID(0)
   {
   }

   void WorkflowBase::checkCores(const int cores)
   {
      if(cores != this->mNcore)
      {
         throw EPMException("WorkflowBase::checkCores", "Workflow and parameters have conflicting number of cores");
      }
   }

}
