/** \file SpectrumFileDefs.hpp
 *  \brief Definition and names used by the spectrum file
 */

#ifndef SPECTRUMFILEDEFS_HPP
#define SPECTRUMFILEDEFS_HPP

// System includes
//
#include <string>

// External includes
//

// Project includes
//

namespace EPMDynamo {

   /**
    * \brief Definition and names used by the spectrum file
    */
   class SpectrumFileDefs
   {
      public:

         /**
          * @brief HEADER part for spectrum file
          */
         static const std::string   HEADER;

         /**
          * @brief VERSION part for spectrum file
          */
         static const std::string   VERSION;

         /**
          * @brief BASENAME of spectrum file
          */
         static const std::string   BASENAME;

         /**
          * @brief EXTENSION of spectrum file
          */
         static const std::string   EXTENSION;

         /**
         * @brief Destructor
         */
         virtual ~SpectrumFileDefs() {};

      private:
         /**
         * @brief Empty destructor
         */
         SpectrumFileDefs();
   };

}

#endif // SPECTRUMFILEDEFS_HPP
