/** \file SerialFFTManipulator.cpp
 *  \brief Source for the serial FFT data manipulator
 */

// System includes
//

// External includes
//

// Class include
//
#include "DataManipulators/SerialFFTManipulator.hpp"

// Project includes
//

namespace EPMDynamo {

   SerialFFTManipulator::SerialFFTManipulator(const SmartTruncation pTrunc)
      : SerialManipulatorBase<FFTFlatScalar, FFTOrderScalar>(pTrunc, 2, 2)
   {
   }

   void SerialFFTManipulator::receive(FFTFlatScalar &rData)
   {
      // Get data from Send queue
      FFTOrderScalar &inMData = this->receiveBData();

      int nM = inMData.nM();
      int nR;
      int nTh = inMData.nTh();

      for(int m = 0; m < nM; ++m)
      {
         nR = inMData.nR(m);
         for(int n = 0; n < nR; ++n)
         {
            for(int th=0; th < nTh; ++th)
            {
               rData.rFlat()(m, n*nTh + th) = inMData.mshell(m)(th,n);
            }
         }
      }

      int padNm = rData.nM();
      for(int m = nM; m < padNm; ++m)
      {
         rData.rFlat().row(m).setConstant(0.0);
      }

      // Free Temporary storage
      freeTmp(inMData);
   }

   void SerialFFTManipulator::receive(FFTOrderScalar &rData)
   {
      // Get data from Send queue
      FFTFlatScalar   &inFData = this->receiveFData();

      int nM = rData.nM();
      int nR;
      int nTh = rData.nTh();

      for(int m = 0; m < nM; ++m)
      {
         nR = rData.nR(m);
         for(int n = 0; n < nR; ++n)
         {
            for(int th=0; th < nTh; ++th)
            {
               rData.rMShell(m)(th,n) = inFData.flat()(m, n*nTh + th) ;
            }
         }
      }

      // Free Temporary storage
      freeTmp(inFData);
   }

}
