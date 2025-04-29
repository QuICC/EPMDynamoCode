/** \file FieldProbeFileDefs.hpp
 *  \brief Definitions and names used by the field probe output files
 */

#ifndef FIELDPROBEFILEDEFS_HPP
#define FIELDPROBEFILEDEFS_HPP

// System includes
//
#include <string>

// External includes
//

// Project includes
//

namespace EPMDynamo {

   /**
    * \brief Definitions and names used by the field probe output files
    */
   class FieldProbeFileDefs
   {
      public:

         /**
          * @brief HEADER part for probe file
          */
         static const std::string   HEADER;

         /**
          * @brief VERSION part for probe file
          */
         static const std::string   VERSION;

         /**
          * @brief BASENAME of probe file
          */
         static const std::string   BASENAME;

         /**
          * @brief EXTENSION of probe file
          */
         static const std::string   EXTENSION;

         /**
         * @brief Destructor
         */
         virtual ~FieldProbeFileDefs(){};

      private:
         /**
         * @brief Empty destructor
         */
         FieldProbeFileDefs();
   };
}

#endif // FIELDPROBEFILEDEFS_HPP
