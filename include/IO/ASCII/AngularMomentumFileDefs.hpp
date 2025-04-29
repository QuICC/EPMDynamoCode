/** \file AngularMomentumFileDefs.hpp
 *  \brief Definitions and names used by the angular momentum diagnostics output files
 */

#ifndef ANGULARMOMENTUMFILEDEFS_HPP
#define ANGULARMOMENTUMFILEDEFS_HPP

// System includes
//
#include <string>

// External includes
//

// Project includes
//

namespace EPMDynamo {

   /**
    * \brief Definitions and names used by the angular momentum diagnostics output files
    */
   class AngularMomentumFileDefs
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
         virtual ~AngularMomentumFileDefs(){};

      private:
         /**
         * @brief Empty destructor
         */
         AngularMomentumFileDefs();
   };
}

#endif // ANGULARMOMENTUMFILEDEFS_HPP
