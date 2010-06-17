/** \file FFTransform.cpp
 *  \brief Source of the implementation of the FFT transform
 */

// System includes
//

// External includes
//

// Class include
//
#include "Transforms/FFT/FFTransform.hpp"

// Project includes
//

namespace EPMDynamo {

   FFTransform::FFTransform(SmartTruncation pTrunc)
      : FFTransformBase(pTrunc)
   {
      this->mTmpRTP.resize(this->mSize, this->mHowmany);

      // Initialise FFTW interface
      this->initFFT();
   }

   FFTransform::~FFTransform()
   {
      // Cleanup memory used by FFTW
      this->cleanupFFT();
   }

   unsigned int FFTransform::mFlags = FFTW_PATIENT;

   void FFTransform::initFFT()
   {
      // Create the two plans
      const int  *fftSize = &this->mSize;

      MatrixZ   tmpCmplx(this->mZRow, this->mHowmany);

      // Create the Real to complex plan
      this->mR2CPlan = fftw_plan_many_dft_r2c(1, fftSize, this->mHowmany, this->mTmpRTP.data(), NULL, 1, this->mSize, reinterpret_cast<fftw_complex* >(tmpCmplx.data()), NULL, 1, this->mZRow, this->mFlags);

      // Create the complex to real plan
      this->mC2RPlan = fftw_plan_many_dft_c2r(1, fftSize, this->mHowmany, reinterpret_cast<fftw_complex* >(tmpCmplx.data()), NULL, 1, this->mZRow, this->mTmpRTP.data(), NULL, 1, this->mSize, this->mFlags);
   }

   void FFTransform::cleanupFFT()
   {
      fftw_destroy_plan(this->mR2CPlan);
      fftw_destroy_plan(this->mC2RPlan);
   }

   void FFTransform::forward(FFTFlatScalar &rFFTVal, const RTPScalar &rtpVal)
   {
      this->flattenRTP(rtpVal);

      // Do transform
      fftw_execute_dft_r2c(this->mR2CPlan, this->mTmpRTP.data(), reinterpret_cast<fftw_complex* >(rFFTVal.rFlat().data()));
   }

   void FFTransform::backward(RTPScalar &rRTPVal, const FFTFlatScalar &fftVal)
   {
      // Do transform
      fftw_execute_dft_c2r(this->mC2RPlan, reinterpret_cast<fftw_complex* >(const_cast<std::complex<double> *>(fftVal.flat().data())), this->mTmpRTP.data());

      this->unflattenRTP(rRTPVal);
   }

   void FFTransform::flattenRTP(const RTPScalar &rtpVal)
   {
      int nR = rtpVal.nR();
      int nPh = rtpVal.nPh();
      int nTh0 = 0;
      int nTh;

      nTh = rtpVal.nTh(0);
      this->mTmpRTP.block(0, 0, nPh, nTh) = this->mScale*rtpVal.shell(0);

      for(int n = 1; n < nR; ++n)
      {
         nTh0 += nTh;
         nTh = rtpVal.nTh(n);
         this->mTmpRTP.block(0, nTh0, nPh, nTh) = this->mScale*rtpVal.shell(n);
      }
   }

   void FFTransform::unflattenRTP(RTPScalar &rRTPVal)
   {
      int nR = rRTPVal.nR();
      int nPh = rRTPVal.nPh();
      int nTh0 = 0;
      int nTh;

      nTh = rRTPVal.nTh(0);
      rRTPVal.rShell(0) = this->mTmpRTP.block(0, 0, nPh, nTh);

      for(int n = 1; n < nR; ++n)
      {
         nTh0 += nTh;
         nTh = rRTPVal.nTh(n);
         rRTPVal.rShell(n) = this->mTmpRTP.block(0, nTh0, nPh, nTh);
      }
   }

}
