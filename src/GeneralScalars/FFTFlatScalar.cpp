/** \file FFTFlatScalar.cpp
 *  \brief Source of the flat storage FFT data type
 */

// System includes
//

// External includes
//

// Class include
//
#include "GeneralScalars/FFTFlatScalar.hpp"

// Project includes
//

namespace EPMDynamo {

   FFTFlatScalar::FFTFlatScalar(const SmartTruncation pTrunc)
      : FFTFlatScalarBase(pTrunc)
   {
      this->initFlatShells();
   }

   void FFTFlatScalar::initFlatShells()
   {
      this->mFlatShells.resize(this->nM(), this->nCols());
   }

}
