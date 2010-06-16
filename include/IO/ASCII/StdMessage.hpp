/** \file StdMessage.hpp
 *  \brief Implementation of std message output
 */

#ifndef STDMESSAGE_HPP
#define STDMESSAGE_HPP

// Configuration includes
//
#include "Config/SmartPointer.h"

// System includes
//
#include <string>

// External includes
//

// Project includes
//
#include "IO/ASCII/ASCIIEWriter.hpp"

namespace EPMDynamo {

   /**
    * @brief Implementation of ASCIIFile for a general message output (replace stdout) ASCII file
    */
   class StdMessage: public ASCIIEWriter
   {
      public:
         /**
         * @brief Constructor 
         *
         * @param name Name of the std message ouput file
         */
         StdMessage(std::string name);

         /**
         * @brief Destructor
         */
         virtual ~StdMessage() {};

         /**
          * @brief Init the file
          */
         void init();

         /**
          * @brief
          */
         virtual void write();

         /**
          * @brief Finalise
          */
         void finalise();
         
      protected:

      private:
         /**
          * Backup std::cout buffer
          */
         std::streambuf   *mpCoutBuffer;
   };

   /// Typedef for a smart pointer of a StdMessage
   typedef EPMSHARED_PTR<StdMessage> SmartStdMessage;

}

#endif // STDMESSAGE_HPP
