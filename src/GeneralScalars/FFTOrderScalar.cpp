/** \file FFTOrderScalar.cpp
 *  \brief Source of the implementation of the harmonic order ordered FFT data type
 */

// System includes
//

// External includes
//

// Class include
//
#include "GeneralScalars/FFTOrderScalar.hpp"

// Project includes
//

namespace EPMDynamo {

   FFTOrderScalar::FFTOrderScalar(const SmartTruncation pTrunc)
      : FFTOrderScalarBase(pTrunc)
   {
      // Init m shells
      initMShells();
   }

   void FFTOrderScalar::initMShells()
   {
      int nM = this->nM();
      for(int m = 0; m < nM; ++m)
      {
         this->mMShells.push_back(MatrixZ(this->nTh(), this->nR(m)));
      }
   }

}
