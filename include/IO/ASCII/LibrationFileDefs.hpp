/** \file LibrationFileDefs.hpp
 *  \brief Definitions and names used by the libration diagnostics output files
 */

#ifndef LIBRATIONFILEDEFS_HPP
#define LIBRATIONFILEDEFS_HPP

// System includes
//
#include <string>

// External includes
//

// Project includes
//

namespace EPMDynamo {

   /**
    * \brief Definitions and names used by the libration diagnostics output files
    */
   class LibrationFileDefs
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
         virtual ~LibrationFileDefs(){};

      private:
         /**
         * @brief Empty destructor
         */
         LibrationFileDefs();
   };
}

#endif // LIBRATIONFILEDEFS_HPP
