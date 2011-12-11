/** \file DipoleFileDefs.hpp
 *  \brief Definitions and names used by the dipole diagnostics output files
 */

#ifndef DIPOLEFILEDEFS_HPP
#define DIPOLEFILEDEFS_HPP

// System includes
//
#include <string>

// External includes
//

// Project includes
//

namespace EPMDynamo {

   /**
    * \brief Definitions and names used by the dipole diagnostics output files
    */
   class DipoleFileDefs
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
         virtual ~DipoleFileDefs(){};

      private:
         /**
         * @brief Empty destructor
         */
         DipoleFileDefs();
   };
}

#endif // DIPOLEFILEDEFS_HPP
