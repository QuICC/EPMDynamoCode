/** \file FDSHShell.cpp
 *  \brief Source of the FDSH spherical shell implementation
 */

// System includes
//

// External includes
//

// Class include
//
#include "GeneralScalars/FDSHShell.hpp"

// Project includes
//

namespace EPMDynamo {

   FDSHShell::FDSHShell(SmartTruncation pTrunc)
      : mpTrunc(pTrunc)
   {
      // Initialise degrees storage
      this->initOrders();
   }

   void FDSHShell::initOrders()
   {
      int nL = this->mpTrunc->sim()->hoz()->nL();
      ArrayI   nLs = this->mpTrunc->local()->fdsh()->nLArray(nL);

      int nM = this->mpTrunc->local()->fdsh()->nM();
      for(int i=0; i < nM; ++i)
      {
         this->mOrders.push_back(SHDegrees(nLs(i)));
      }
   }

}
