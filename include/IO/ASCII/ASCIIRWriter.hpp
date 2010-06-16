/** \file ASCIIRWriter.hpp
 *  \brief Implementation of replacing ASCII writer
 */

#ifndef ASCIIRWRITER_HPP
#define ASCIIRWRITER_HPP

// System includes
//
#include <fstream>

// External includes
//

// Project includes
//
#include "IO/ASCII/ASCIIWriter.hpp"

namespace EPMDynamo {

   /**
    * \brief Implementation of replacing ASCII writer
    */
   class ASCIIRWriter: public ASCIIWriter
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
         ASCIIRWriter(std::string name, std::string ext, std::string header, std::string version);

         /**
         * @brief Destructor
         */
         virtual ~ASCIIRWriter() {};

         /**
          * @brief Initialise the file
          */
         virtual void init();

         /**
          * @brief Write the content
          */
         virtual void write() = 0;

         /**
          * @brief Finalise the file
          */
         virtual void finalise();
         
      protected:

         /**
          * @brief Operation to perform just before writing data
          */
         void preWrite();

         /**
          * @brief Operation to perform just after writing data
          */
         void postWrite();

      private:
   };
}

#endif // ASCIIRWRITER_HPP
