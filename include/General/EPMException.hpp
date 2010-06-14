/** \file EPMException.hpp
 *  \brief Simple exception class.
 */

#ifndef EPMEXCEPTION_HPP
#define EPMEXCEPTION_HPP

// System includes
//
#include <string>

// Project includes
//

namespace EPMDynamo {

   /**
    * @brief Simple exception class
    *
    * This is very simple exception class. It allows to give a location and a message string.
    * These can then be recovered in a catch block.
    */
   class EPMException
   {
      public:
         /**
         * @brief Constructs the exception based on the given strings.
         *
         * \param location Location where the exception has been thrown 
         * \param message Error message
         */
         EPMException(std::string location = "Everywhere", std::string message = "Aaaarghh!");

         /**
         * @brief Simple empty destructor
         */
         virtual ~EPMException() {};

         /**
          * @brief Get the location message
          */
         std::string location() const;

         /**
          * @brief Get the error message
          */
         std::string error() const;

         /**
          * @brief Print standard error message
          *
          * Prints the location and error message in a nice for to std::cerr stream.
          */
         void printStdMessage() const;
         
      protected:

      private:
         /**
          * @brief Storage for the location string
          */
         std::string    mLocation;

         /**
          * @brief Storage for the error message string
          */
         std::string    mMessage;
   };

   inline std::string EPMException::location() const
   {
      return this->mLocation;
   }

   inline std::string EPMException::error() const
   {
      return this->mMessage;
   }

}

#endif // EPMEXCEPTION_HPP
