/** \file FFTOrderScalar.hpp
 *  \brief Implementation of the FFT harmonic order ordered data type
 */

#ifndef FFTORDERSCALAR_HPP
#define FFTORDERSCALAR_HPP

// System includes
//
#include <vector>

// External includes
//

// Project includes
//
#include "General/EPMTypedefs.hpp"
#include "Domain/Truncation.hpp"
#include "GeneralScalars/FFTOrderScalarBase.hpp"

namespace EPMDynamo {

   /**
    * \brief Implementation of the FFT harmonic order ordered data type
    */
   class FFTOrderScalar: public FFTOrderScalarBase
   {
      public:
         /**
          * @brief Constructor
          *
          * @param pTrunc Truncation information pointer
          */
         FFTOrderScalar(const SmartTruncation pTrunc);

         /**
          * @brief Destructor
          */
         virtual ~FFTOrderScalar() {};

         /**
          * @brief Get the harmonic order m "Harmonic shell"
          *
          * @param m Harmonic order m
          */
         const MatrixZ&    mshell(const int m) const;

         /**
          * @brief Set the harmonic order m "Harmonic shell"
          *
          * @param m Harmonic order m
          */
         MatrixZ&    rMShell(const int m);
         
      protected:

      private:
         /**
          * @brief Storage for the values in MShell ordering
          */
         std::vector<MatrixZ>  mMShells;

         /**
          * @brief Initialise the mshells
          */
         void initMShells();
   };

   inline const MatrixZ&   FFTOrderScalar::mshell(const int m) const
   {
      return this->mMShells.at(m);
   }

   inline MatrixZ&   FFTOrderScalar::rMShell(const int m)
   {
      return this->mMShells.at(m);
   }

}

#endif // FFTORDERSCALAR_HPP
