/** \file XMLReaderDefs.hpp
 *  \brief Definitions and names used by the XML reader
 */

#ifndef XMLREADERDEFS_HPP
#define XMLREADERDEFS_HPP

// System includes
//
#include <fstream>

// External includes
//

// Project includes
//
#include "IO/ASCII/XMLReaderDefs.hpp"

namespace EPMDynamo {

   /**
    * @brief Simpler reader class for XML data
    *
    * This will mainly be used for reading in the parameters
    */
   class XMLReaderDefs
   {
      public:
         /**
          * @brief HEADER part for parameters file
          */
         static const std::string   GENERALXML;

         /**
          * @brief VERSION part for parameters file
          */
         static const std::string   HEADERXML;

         /**
          * @brief BASENAME of parameters file
          */
         static const std::string   VERSIONXML;

         /**
         * @brief Destructor
         */
         virtual ~XMLReaderDefs() {};
         
      protected:
         /**
         * @brief Constructor
         */
         XMLReaderDefs();

      private:
   };

}

#endif // XMLREADERDEFS_HPP
