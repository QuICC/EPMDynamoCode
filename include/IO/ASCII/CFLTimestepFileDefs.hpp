/** \file CFLTimestepFileDefs.hpp 
 *  \brief Definition and names used by the time information file
 */

#ifndef CFLTIMESTEPFILEDEFS_HPP
#define CFLTIMESTEPFILEDEFS_HPP

// System includes
//
#include <string>

// External includes
//

// Project includes
//

namespace EPMDynamo {

   /**
    * \brief Definition and names used by the time information file
    */
   class CFLTimestepFileDefs
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
         virtual ~CFLTimestepFileDefs() {};

      private:
         /**
         * @brief Empty destructor
         */
         CFLTimestepFileDefs();
   };

}

#endif // CFLTIMESTEPFILEDEFS_HPP
