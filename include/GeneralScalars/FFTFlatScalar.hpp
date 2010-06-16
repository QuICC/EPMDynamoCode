/** \file FFTFlatScalar.hpp
 *  \brief Implementation of the FFT flat storage data type
 */

#ifndef FFTFLATSCALAR_HPP
#define FFTFLATSCALAR_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "General/EPMTypedefs.hpp"
#include "Domain/Truncation.hpp"
#include "GeneralScalars/FFTFlatScalarBase.hpp"

namespace EPMDynamo {

   /**
    * \brief Implementation of the FFT flat storage data type
    */
   class FFTFlatScalar: public FFTFlatScalarBase
   {
      public:
         /**
          * @brief Constructor
          *
          * @param pTrunc Truncation information
          */
         FFTFlatScalar(const SmartTruncation pTrunc);

         /**
          * @brief Destructor
          */
         virtual ~FFTFlatScalar() {};

         /**
          * @brief No docu
          *
          * \todo Need documentation
          */
         const MatrixZ&    flat() const;

         /**
          * @brief No docu
          *
          * \todo Need documentation
          */
         MatrixZ&    rFlat();

         /**
          * @brief Zero padding for FFT
          */
         void doZeroPadding();
         
      protected:

      private:
         /**
          * @brief Storage for the flat values
          */
         MatrixZ  mFlatShells;

         /**
          * @brief Initialise the flat shells
          */
         void initFlatShells();
   };

   inline const MatrixZ&   FFTFlatScalar::flat() const
   {
      return this->mFlatShells;
   }

   inline MatrixZ&   FFTFlatScalar::rFlat()
   {
      return this->mFlatShells;
   }

   inline void FFTFlatScalar::doZeroPadding()
   {
      int nPad = this->nM() - this->mpTrunc->sim()->hoz()->nM();
      this->mFlatShells.block(this->mpTrunc->sim()->hoz()->nM(), 0, nPad, this->nCols()).setConstant(0.0);
   }

}

#endif // FFTFLATSCALAR_HPP
