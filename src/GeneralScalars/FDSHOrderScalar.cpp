/** \file FDSHOrderScalar.cpp
 *  \brief Source of the implementation of the harmonic order ordered FDSH data type
 */

// System includes
//

// External includes
//

// Class include
//
#include "GeneralScalars/FDSHOrderScalar.hpp"

// Project includes
//

namespace EPMDynamo {

   FDSHOrderScalar::FDSHOrderScalar(const SmartTruncation pTrunc)
      : FDSHOrderScalarBase(pTrunc)
   {
      // Init the SH degrees
      this->initMShells();
   }

   void FDSHOrderScalar::initMShells()
   {
      int nM = this->nM();
      for(int m = 0; m < nM; ++m)
      {
         this->mMShells.push_back(MatrixZ(this->nL(m), this->nR(m)));
      }
   }

}
