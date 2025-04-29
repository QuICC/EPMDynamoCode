/** \file ETDInfluenceKernel.cpp
 *  \brief Implemenation of influence of the kernel solution
 */


// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "Timestepping/ETD/ETDInfluenceKernel.hpp"

// Project includes
//

namespace EPMDynamo {

   ETDInfluenceKernel::ETDInfluenceKernel(EPMFloat nFactor, SmartETDOperators pOpM1)
      : mNFactor(nFactor), mpOpM1(pOpM1)
   {
   }

   void ETDInfluenceKernel::computeInfluence(Array &rKernel, const int l)
   {
      int rows = this->mpOpM1->harmOp(l).op().rows();
      rKernel.topRows(rows) = this->mNFactor*this->mpOpM1->harmOp(l).op() * rKernel.topRows(rows);
      this->mpOpM1->harmOp(l).extend(rKernel);
   }

}
