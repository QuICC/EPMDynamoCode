/** \file ASCIIFile.hpp
 *  \brief Implementation of a general ASCIIFile
 */

#ifndef ASCIIFILE_HPP
#define ASCIIFILE_HPP

// System includes
//
#include <string>
#include <vector>

// External includes
//

// Project includes
//
#include "IO/FileBase.hpp"

namespace EPMDynamo {

   /**
    * @brief Implementation of FileBase for a general ASCII file
    */
   class ASCIIFile: public FileBase
   {
      public:
         /**
         * @brief Constructor
         *
         * @param name Filename
         * @param ext File extension
         * @param header Header string of file
         * @param version  Version string of file 
         */
         ASCIIFile(std::string name, std::string ext, std::string header, std::string version);

         /**
         * @brief Destructor
         */
         virtual ~ASCIIFile() {};
         
      protected:

         /**
          * @brief Header of file before header information
          */
         static const std::string FILE_HEADER;

         /**
          * @brief Header of file before version information
          */
         static const std::string FILE_VERSION;

         /**
          * @brief Header of the file to check compatibility
          */
         std::string mHeader;

         /**
          * @brief Version of the file to check compatibility
          */
         std::string mVersion;

      private:
   };

}

#endif // ASCIIFILE_HPP
