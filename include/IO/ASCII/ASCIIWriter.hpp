/** \file ASCIIWriter.hpp
 *  \brief Implementation of a general ASCII writer
 */

#ifndef ASCIIWRITER_HPP
#define ASCIIWRITER_HPP

// Configuration includes
//
#include "Config/SmartPointer.h"

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
    * \brief Implementation of a general ASCII writer
    */
   class ASCIIWriter: public ASCIIFile
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
         ASCIIWriter(std::string name, std::string ext, std::string header, std::string version);

         /**
         * @brief Destructor
         */
         virtual ~ASCIIWriter() {};

         /**
          * @brief Initialise the file
          */
         virtual void init() = 0;

         /**
          * @brief Write the content
          */
         virtual void write() = 0;

         /**
          * @brief Finalise the file
          */
         virtual void finalise() = 0;
         
      protected:
         /**
          * @brief Handle to the file
          */
         std::ofstream mFile;

         /**
          * @brief Open the file
          */
         void open();

         /**
          * @brief Close the file
          */
         void close();

      private:
   };

   /// Typedef for a smart pointer of a ASCIIWriter
   typedef EPMSHARED_PTR<ASCIIWriter> SmartASCIIWriter;
}

#endif // ASCIIWRITER_HPP
