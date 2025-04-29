/** \file ASCIINWriter.hpp
 *  \brief Implementation of a numbered ASCII file writer
 */

#ifndef ASCIINWRITER_HPP
#define ASCIINWRITER_HPP

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
    * @brief Implementation of a numbered ASCII file writer
    */
   class ASCIINWriter: public ASCIIWriter
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
         ASCIINWriter(std::string name, std::string ext, std::string header, std::string version);

         /**
         * @brief Destructor
         */
         virtual ~ASCIINWriter() {};

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

         /**
          * @brief Update the file name with counter value
          */
         void updateName();

      private:
         /**
          * @brief Zero Fill width
          */
         static const int msIDWidth;

         /**
          * @brief File counter
          */
         int mCounter;

         /**
          * @brief Base name of the file
          */
         const std::string    mBaseName;
   };
}

#endif // ASCIINWRITER_HPP
