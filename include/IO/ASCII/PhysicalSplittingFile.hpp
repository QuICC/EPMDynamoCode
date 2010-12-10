/** \file PhysicalSplittingFile.hpp
 *  \brief Implementation of an ASCII physical space splitting description file
 */

#ifndef PHYSICALSPLITTINGFILE_HPP
#define PHYSICALSPLITTINGFILE_HPP

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
    * \brief Implementation of an ASCII physical space splitting description file
    */
   class PhysicalSplittingFile: public SplittingFile
   {
      public:
         /**
         * @brief Constructor
         *
         * @param pTrunc Smart truncation information
         */
         PhysicalSplittingFile(const SmartTruncation pTrunc);

         /**
         * @brief Destructor
         */
         virtual ~PhysicalSplittingFile() {};

         /**
          * @brief Write the physical space splitting information to file
          */
         virtual void write();
         
      protected:

      private:
   };
}

#endif // PHYSICALSPLITTINGFILE_HPP
