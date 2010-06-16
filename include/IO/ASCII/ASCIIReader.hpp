/** \file ASCIIReader.hpp
 *  \brief Implementation of an ASCII file reader
 */

#ifndef ASCIIREADER_HPP
#define ASCIIREADER_HPP

// System includes
//
#include <fstream>

// External includes
//

// Project includes
//
#include "IO/ASCII/ASCIIFile.hpp"

namespace EPMDynamo {

   /**
    * @brief Implementation of basic ASCII file reader
    *
    * This will mainly be used for reading in the parameters
    */
   class ASCIIReader: public ASCIIFile
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
         ASCIIReader(std::string name, std::string ext, std::string header, std::string version);

         /**
         * @brief Destructor
         */
         virtual ~ASCIIReader() {};

         /**
          * @brief Initialise the file
          */
         virtual void init();

         /**
          * @brief Read the content
          */

         virtual void read() = 0;

         /**
          * @brief Finalise the file
          */
         virtual void finalise();
         
      protected:
         /**
          * @brief Handle to the file
          */
         std::ifstream mFile;

         /**
          * @brief Open the file
          */
         void open();

         /**
          * @brief Close the file
          */
         void close();

         /**
          * @brief Check compatibility of opened file
          */
         virtual void checkCompatibility();

      private:
   };

}

#endif // ASCIIREADER_HPP
