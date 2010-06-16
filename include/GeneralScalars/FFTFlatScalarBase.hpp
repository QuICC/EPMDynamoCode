/** \file FFTFlatScalarBase.hpp
 *  \brief Base of the FFT flat storage data type implementation
 */

#ifndef FFTFLATSCALARBASE_HPP
#define FFTFLATSCALARBASE_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "General/EPMTypedefs.hpp"
#include "Domain/Truncation.hpp"

namespace EPMDynamo {

   /**
    * \brief Base of the FFT flat storage data type implementation
    */
   class FFTFlatScalarBase
   {
      public:
         /**
          * @brief Constructor
          *
          * @param pTrunc Truncation information
          */
         FFTFlatScalarBase(const SmartTruncation pTrunc);

         /**
          * @brief Destructor
          */
         virtual ~FFTFlatScalarBase() {};

         /**
          * @brief Get the number of harmonic orders 
          */
         int   nM() const;

         /**
          * @brief Get the number of theta points for given radial index 
          *
          * @param n Radial index n
          */
         int   nTh(const int n) const;

         /**
          * @brief Get the number of radii 
          */
         int   nR() const;

         /**
          * @brief Get total number of columns 
          */
         int   nCols() const;
         
      protected:
         /**
          * @brief Spectral truncation information
          */
         SmartTruncation   mpTrunc;

      private:
   };

   inline int FFTFlatScalarBase::nM() const
   {
      return this->mpTrunc->sim()->hoz()->nPh()/2+1;
   }

   inline int   FFTFlatScalarBase::nR() const
   {
      return this->mpTrunc->local()->rtp()->nR();
   }

   inline int   FFTFlatScalarBase::nTh(const int n) const
   {
      return this->mpTrunc->local()->rtp()->nTh(n);
   }

   inline int   FFTFlatScalarBase::nCols() const
   {
      return this->mpTrunc->local()->rtp()->nThArray().sum();
   }

}

#endif // FFTFLATSCALARBASE_HPP
