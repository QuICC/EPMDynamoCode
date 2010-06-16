/** \file SourceDefs.hpp
 *  \brief Definitions and names used by the source field readers/writers
 */

#ifndef SOURCEDEFS_HPP
#define SOURCEDEFS_HPP

// System includes
//
#include <string>

// External includes
//

// Project includes
//

namespace EPMDynamo {

   /**
    * \brief Definitions and names used by the source field readers/writers
    */
   class SourceDefs
   {
      public:
         /**
          * @brief HEADER part for State file
          */
         static const std::string   TOROIDALTAG;

         /**
          * @brief HEADER part for State file
          */
         static const std::string   POLOIDALTAG;

         /**
          * @brief HEADER part for State file
          */
         static const std::string   CODENSITYTAG;

         /**
          * @brief HEADER part for State file
          */
         static const std::string   MAGNETICTAG;

         /**
          * @brief HEADER part for State file
          */
         static const std::string   VELOCITYTAG;

         /**
          * @brief HEADER part for State file
          */
         static const std::string   HEADER;

         /**
          * @brief VERSION part for State file
          */
         static const std::string   VERSION;

         /**
          * @brief BASENAME of State file
          */
         static const std::string   BASENAME;

         /**
          * @brief EXTENSION of State file
          */
         static const std::string   EXTENSION;

         /**
         * @brief Destructor
         */
         virtual ~SourceDefs() {};

      private:
         /**
         * @brief Empty destructor
         */
         SourceDefs();
   };
}

#endif // SOURCEDEFS_HPP
