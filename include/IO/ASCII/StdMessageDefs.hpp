/** \file StdMessageDefs.hpp
 *  \brief Definitions and names used by the Std message output
 */

#ifndef STDMESSAGEDEFS_HPP
#define STDMESSAGEDEFS_HPP

// System includes
//
#include <string>

// External includes
//

// Project includes
//

namespace EPMDynamo {

   /**
    * \brief Definitions and names used by the Std message output
    */
   class StdMessageDefs
   {
      public:

         /**
          * @brief HEADER part for StdMessage file
          */
         static const std::string   HEADER;

         /**
          * @brief VERSION part for StdMessage file
          */
         static const std::string   VERSION;

         /**
          * @brief BASENAME of StdMessage file
          */
         static const std::string   BASENAME;

         /**
          * @brief EXTENSION of StdMessage file
          */
         static const std::string   EXTENSION;

         /**
         * @brief Destructor
         */
         virtual ~StdMessageDefs() {};

      private:
         /**
         * @brief Empty destructor
         */
         StdMessageDefs();
   };

}

#endif // STDMESSAGEDEFS_HPP
