/** \file ConfigurationFileDefs.hpp
 *  \brief Definitions and names used by the configuration file
 */

#ifndef CONFIGURATIONFILEDEFS_HPP
#define CONFIGURATIONFILEDEFS_HPP

// System includes
//
#include <string>

// External includes
//

// Project includes
//

namespace EPMDynamo {

   /**
    * \brief Definitions and names used by the configuration file
    */
   class ConfigurationFileDefs
   {
      public:

         /**
          * @brief HEADER part for parameters file
          */
         static const std::string   HEADER;

         /**
          * @brief VERSION part for parameters file
          */
         static const std::string   VERSION;

         /**
          * @brief BASENAME of parameters file
          */
         static const std::string   BASENAME;

         /**
          * @brief EXTENSION of parameters file
          */
         static const std::string   EXTENSION;

         /**
          * @brief Truncation part of parameters file
          */
         static const std::string   CONFIGXML;

         /**
         * @brief Destructor
         */
         virtual ~ConfigurationFileDefs() {};

      private:
         /**
         * @brief Empty destructor
         */
         ConfigurationFileDefs();
   };
}

#endif // CONFIGURATIONFILEDEFS_HPP
