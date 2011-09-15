/** \file ExtractRadiusFileDefs.hpp
 *  \brief Definitions and names used by the extract radius diagnostics output files
 */

#ifndef EXTRACTRADIUSFILEDEFS_HPP
#define EXTRACTRADIUSFILEDEFS_HPP

// System includes
//
#include <string>

// External includes
//

// Project includes
//

namespace EPMDynamo {

   /**
    * \brief Definitions and names used by the extract radius diagnostics output files
    */
   class ExtractRadiusFileDefs
   {
      public:

         /**
          * @brief HEADER part for energy file
          */
         static const std::string   HEADER;

         /**
          * @brief VERSION part for energy file
          */
         static const std::string   VERSION;

         /**
          * @brief BASENAME of energy file
          */
         static const std::string   BASENAME;

         /**
          * @brief EXTENSION of energy file
          */
         static const std::string   EXTENSION;

         /**
         * @brief Destructor
         */
         virtual ~ExtractRadiusFileDefs(){};

      private:
         /**
         * @brief Empty destructor
         */
         ExtractRadiusFileDefs();
   };
}

#endif // EXTRACTRADIUSFILEDEFS_HPP
