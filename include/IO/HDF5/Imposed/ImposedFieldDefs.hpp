/** \file ImposedFieldDefs.hpp
 *  \brief Definitions and names used by the imposed field readers/writers
 */

#ifndef IMPOSEDFIELDDEFS_HPP
#define IMPOSEDFIELDDEFS_HPP

// System includes
//
#include <string>

// External includes
//

// Project includes
//

namespace EPMDynamo {

   /**
    * \brief Definitions and names used by the imposed field readers/writers
    */
   class ImposedFieldDefs
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
         virtual ~ImposedFieldDefs() {};

      private:
         /**
         * @brief Empty destructor
         */
         ImposedFieldDefs();
   };
}

#endif // IMPOSEDFIELDDEFS_HPP
