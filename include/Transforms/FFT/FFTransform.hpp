/** \file FFTransform.hpp
 *  \brief Implementation of the FFT transform
 */

#ifndef FFTRANSFORM_HPP
#define FFTRANSFORM_HPP

// System includes
//
#include <complex>
#include <fftw3.h>

// External includes
//

// Project includes
//
#include "General/EPMTypedefs.hpp"
#include "GeneralScalars/RTPScalar.hpp"
#include "GeneralScalars/FFTFlatScalar.hpp"
#include "Transforms/FFT/FFTransformBase.hpp"

namespace EPMDynamo {

   /**
    * @brief Interface class to the FFTW routines
    *
    * This class is reponsible for performing the FFT. It is written
    * as an interface to the FFTW functions. The transforms a performed with the help
    * of two variable inside the class. This will require care on a later step with
    * MPI/threading.
    */
   class FFTransform : public FFTransformBase
   {
      public:
         /**
         * @brief Constructs the required matrices and plans for the FFT computation
         *
         * This constructor will initialise the required array and set up the
         * forward and backward transform plans.
         *
         * @param pTrunc Truncation information
         */
         FFTransform(SmartTruncation pTrunc);

         /**
         * @brief The FFTW plans will be destroyed here
         */
         virtual ~FFTransform();

         /**
          * @brief Do FFT from real space to complex
          *
          * @param rFFTVal Output FFT transformed values
          * @param rtpVal Input RTP values
          */
         void forward(FFTFlatScalar& rFFTVal, const RTPScalar &rtpVal);

         /**
          * @brief Do FFT from complex space to real
          *
          * @brief rRTPVal Output RTP values
          * @brief fftVal Input FFT values
          */
         void backward(RTPScalar &rRTPVal, const FFTFlatScalar& fftVal);
         
      protected:

      private:
         /**
          * @brief FFTW3 flags for the transforms
          *
          * This is currently simply set to perform a thorough analysis of the 
          * algorithm to use
          */
         static unsigned int  mFlags;

         /**
          * @brief FFTW plan for the real to complex transform (forward)
          */
         fftw_plan   mR2CPlan;

         /**
          * @brief FFTW plan for the complex to real transform (backward)
          */
         fftw_plan   mC2RPlan;

         /**
          * @brief Flat RTP data storage
          */
         Matrix   mTmpRTP;

         /**
          * @brief Initialise the FFTW transforms (i.e. create plans, etc)
          */
         void initFFT();

         /**
          * @brief Cleanup memory used by FFTW on destruction
          */
         void cleanupFFT();

         /**
          * @brief Flatten RTP data
          *
          * @param rtpVal RTP values to flatten
          */
         void flattenRTP(const RTPScalar &rtpVal);

         /**
          * @brief Unflatten RTP data
          *
          * @param rRTPVal Storage for the unflattened values
          */
         void unflattenRTP(RTPScalar &rRTPVal);
   };

}

#endif // FFTRANSFORM_HPP
