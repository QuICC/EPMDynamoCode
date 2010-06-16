/** \file FFTOrderScalarBase.hpp
 *  \brief Base of the FFT harmonic order ordered data type implementation
 */

#ifndef FFTORDERSCALARBASE_HPP
#define FFTORDERSCALARBASE_HPP

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
    * \brief Base of the FFT harmonic order ordered data type implementation
    */
   class FFTOrderScalarBase
   {
      public:
         /**
          * @brief Constructor
          *
          * @param pTrunc Truncation information
          */
         FFTOrderScalarBase(const SmartTruncation pTrunc);

         /**
          * @brief Destructor
          */
         virtual ~FFTOrderScalarBase() {};

         /**
          * @brief Get the number of theta points 
          */
         int   nTh() const;

         /**
          * @brief Get the number of harmonic orders 
          */
         int   nM() const;

         /**
          * @brief Get the number of radial shells 
          *
          * @param m Harmonic order m
          */
         int   nR(const int m) const;

         /**
          * @brief Get the number of radial shells 
          */
         const ArrayI&   nR() const;
         
      protected:
         /**
          * @brief Spectral truncation information
          */
         SmartTruncation   mpTrunc;

      private:
   };

   inline int FFTOrderScalarBase::nTh() const
   {
      return this->mpTrunc->sim()->hoz()->nTh();
   }

   inline int FFTOrderScalarBase::nM() const
   {
      return this->mpTrunc->local()->fdsh()->nM();
   }

   inline int FFTOrderScalarBase::nR(const int m) const
   {
      return this->mpTrunc->local()->fdsh()->nR(m);
   }

   inline const ArrayI& FFTOrderScalarBase::nR() const
   {
      return this->mpTrunc->local()->fdsh()->nRArray();
   }

}

#endif // FFTORDERSCALARBASE_HPP
