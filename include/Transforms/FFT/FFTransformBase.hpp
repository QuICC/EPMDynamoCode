/** \file FFTransformBase.hpp
 *  \brief Base of the implementation of the FFT transforms
 */

#ifndef FFTRANSFORMBASE_HPP
#define FFTRANSFORMBASE_HPP

// System includes
//
#include <complex>

// External includes
//

// Project includes
//
#include "General/EPMTypedefs.hpp"
#include "Domain/Truncation.hpp"
#include "GeneralScalars/RTPFlatScalar.hpp"

namespace EPMDynamo {

   /**
    * \brief Base of the implementation of the FFT transforms
    */
   class FFTransformBase
   {
      public:
         /**
         * @brief Constructor
         *
         * @param pTrunc Truncation information
         */
         FFTransformBase(SmartTruncation pTrunc);

         /**
         * @brief Destructor
         */
         virtual ~FFTransformBase() {};
         
      protected:
         /**
          * @brief Basis size of the transform
          */
         int mSize;

         /**
          * @brief Length of a row of the complex array
          */
         int mZRow;

         /**
          * @brief Number of transforms to perform
          */
         int mHowmany;

         /**
          * @brief Size of the real array for the transform
          */
         int mReSize;

         /**
          * @brief Size of the complex array for the transform
          */
         int mZSize;

         /**
          * @brief Scaling factor for the transform
          *
          * FFTW performs unscale transforms and have thus to be rescaled.
          */
         double mScale;

      private:
         /**
          * @brief Initialise the different sizes
          */
         void initSizes();
   };

}

#endif // FFTRANSFORMBASE_HPP
