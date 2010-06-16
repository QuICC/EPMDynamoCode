/** \file FDSHDegreeScalar.cpp
 *  \brief Source of the implementation of the harmonic degree ordered data type
 */

// System includes
//

// External includes
//

// Class include
//
#include "GeneralScalars/FDSHDegreeScalar.hpp"

// Project includes
//

namespace EPMDynamo {

   FDSHDegreeScalar::FDSHDegreeScalar(const SmartTruncation pTrunc)
      : FDSHDegreeScalarBase(pTrunc)
   {
      // Init the SH orders
      this->initLShells();
   }

   void FDSHDegreeScalar::initLShells()
   {
      int nL = this->nL();
      for(int l = 0; l < nL; ++l)
      {
         mLShells.push_back(MatrixZ(this->nR(), this->nM(l)));
      }
   }

}
