/** \file SerialWorkflow.cpp
 *  \brief Source of implementation of serial workflow control
 */

// System includes
//

// External includes
//

// Class include
//
#include "Parallelisers/SerialWorkflow.hpp"

// Project includes
//

namespace EPMDynamo {

   SerialWorkflow::SerialWorkflow(const int nCore)
      : WorkflowBase(nCore)
   {
      this->setup();
   }

   void SerialWorkflow::init()
   {
   }

   void SerialWorkflow::setup()
   {
      this->checkCores(1);
   }

   void SerialWorkflow::finalize()
   {
   }

}
