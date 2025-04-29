/** \file TimeFileDefs.hpp 
 *  \brief Definition and names used by the time information file
 */

#ifndef TIMEFILEDEFS_HPP
#define TIMEFILEDEFS_HPP

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
   class TimeFileDefs
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
         virtual ~TimeFileDefs() {};

      private:
         /**
         * @brief Empty destructor
         */
         TimeFileDefs();
   };

}

#endif // TIMEFILEDEFS_HPP
