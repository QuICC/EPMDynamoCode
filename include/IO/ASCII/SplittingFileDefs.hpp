/** \file SplittingFileDefs.hpp
 *  \brief Definition and names used by the splitting file
 */

#ifndef SPLITTINGFILEDEFS_HPP
#define SPLITTINGFILEDEFS_HPP

// System includes
//
#include <string>

// External includes
//

// Project includes
//

namespace EPMDynamo {

   /**
    * \brief Definition and names used by the splitting file
    */
   class SplittingFileDefs
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
         virtual ~SplittingFileDefs() {};

      private:
         /**
         * @brief Empty destructor
         */
         SplittingFileDefs();
   };

}

#endif // SPLITTINGFILEDEFS_HPP
