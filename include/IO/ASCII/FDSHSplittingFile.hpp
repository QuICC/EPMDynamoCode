/** \file FDSHSplittingFile.hpp
 *  \brief Implementation of an ASCII FDSH space splitting description file
 */

#ifndef FDSHSPLITTINGFILE_HPP
#define FDSHSPLITTINGFILE_HPP

// System includes
//
#include <string>

// External includes
//

// Project includes
//
#include "Domain/Truncation.hpp"
#include "IO/ASCII/SplittingFile.hpp"

namespace EPMDynamo {

   /**
    * \brief Implementation of an ASCII FDSH space splitting description file
    */
   class FDSHSplittingFile: public SplittingFile
   {
      public:
         /**
         * @brief Constructor
         *
         * @param pTrunc Smart truncation information
         */
         FDSHSplittingFile(const SmartTruncation pTrunc);

         /**
         * @brief Destructor
         */
         virtual ~FDSHSplittingFile() {};

         /**
          * @brief Write the physical space splitting information to file
          */
         virtual void write();
         
      protected:

      private:
   };
}

#endif // FDSHSPLITTINGFILE_HPP
