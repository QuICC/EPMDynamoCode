/** \file SpectralFileDefs.hpp
 *  \brief Definitions and names used by the spectral space data files
 */

#ifndef SPECTRALFILEDEFS_HPP
#define SPECTRALFILEDEFS_HPP

// System includes
//
#include <string>

// External includes
//

// Project includes
//

namespace EPMDynamo {

   /**
    * @brief Implementation of HDF5File for a spectral data HDF5 file
    */
   class SpectralFileDefs
   {
      public:
         /**
          * @brief HEADER part for State file
          */
         static const std::string   TRUNCATIONTAG;

         /**
          * @brief HEADER part for State file
          */
         static const std::string   TRUNCNTAG;

         /**
          * @brief HEADER part for State file
          */
         static const std::string   TRUNCLTAG;

         /**
          * @brief HEADER part for State file
          */
         static const std::string   TRUNCMTAG;

         /**
          * @brief HEADER part for State file
          */
         static const std::string   TRUNCMPTAG;

         /**
         * @brief Destructor
         */
         virtual ~SpectralFileDefs() {};

      private:
         /**
         * @brief Empty destructor
         */
         SpectralFileDefs();
   };
}

#endif // SPECTRALFILEDEFS_HPP
