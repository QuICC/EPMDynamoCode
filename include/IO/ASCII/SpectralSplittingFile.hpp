/** \file SpectralSplittingFile.hpp
 *  \brief Implementation of an ASCII spectral space splitting description file
 */

#ifndef SPECTRALSPLITTINGFILE_HPP
#define SPECTRALSPLITTINGFILE_HPP

// System includes
//
#include <string>

// External includes
//

// Project includes
//
#include "Domain/Truncation.hpp"
#include "IO/ASCII/SplittingFile.hpp"

namespace EPMDynamo {

   /**
    * \brief Implementation of an ASCII spectral space splitting description file
    */
   class SpectralSplittingFile: public SplittingFile
   {
      public:
         /**
         * @brief Constructor
         *
         * @param pTrunc Smart truncation information
         */
         SpectralSplittingFile(const SmartTruncation pTrunc);

         /**
         * @brief Destructor
         */
         virtual ~SpectralSplittingFile() {};

         /**
          * @brief Write the physical space splitting information to file
          */
         virtual void write();
         
      protected:

      private:
   };
}

#endif // SPECTRALSPLITTINGFILE_HPP
