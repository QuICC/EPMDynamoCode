/** \file SplittingFile.hpp
 *  \brief Implementation of an ASCII splitting description file
 */

#ifndef SPLITTINGFILE_HPP
#define SPLITTINGFILE_HPP

// System includes
//
#include <string>

// External includes
//

// Project includes
//
#include "Domain/Truncation.hpp"
#include "IO/ASCII/SplittingFileDefs.hpp"
#include "IO/ASCII/ASCIIRWriter.hpp"

namespace EPMDynamo {

   /**
    * \brief Implementation of an ASCII splitting description file
    */
   class SplittingFile: public ASCIIRWriter
   {
      public:
         /**
         * @brief Constructor
         *
         * @param pTrunc Smart truncation information
         * @param name File name
         */
         SplittingFile(const SmartTruncation pTrunc, std::string name);

         /**
         * @brief Destructor
         */
         virtual ~SplittingFile() {};

         /**
          * @brief Write the splitting information to file
          */
         virtual void write() = 0;
         
      protected:

         /**
          * @brief Storage for the smart truncation information
          */
         SmartTruncation mpTrunc;

      private:
   };
}

#endif // SPLITTINGFILE_HPP
