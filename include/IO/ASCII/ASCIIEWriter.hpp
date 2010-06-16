/** \file ASCIIEWriter.hpp
 *  \brief General implementation of an ASCII "extending" writer
 */

#ifndef ASCIIEWRITER_HPP
#define ASCIIEWRITER_HPP

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
    * @brief Small ASCII file writing class
    *
    * It will be used to output some useful diagnostics values during runtime
    */
   class ASCIIEWriter: public ASCIIWriter
   {
      public:
         /**
         * @brief Constructor
         *
         * @param name Filename
         * @param ext File extension
         * @param header Header string of file
         * @param version Version string of file
         */
         ASCIIEWriter(std::string name, std::string ext, std::string header, std::string version);

         /**
         * @brief Destructor
         */
         virtual ~ASCIIEWriter() {};

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

#endif // ASCIIEWRITER_HPP
