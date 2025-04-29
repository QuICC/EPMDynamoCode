/** \file ConfigurationFile.hpp 
 *  \brief Implementation of a configuration file reader
 */

#ifndef CONFIGURATIONFILE_HPP
#define CONFIGURATIONFILE_HPP

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
#include "General/EPMTypedefs.hpp"
#include "IO/ASCII/XMLReader.hpp"

namespace EPMDynamo {

   /**
    * \brief Implementation of a configuration file reader
    */
   class ConfigurationFile: public XMLReader
   {
      public:
         /**
         * @brief Constructor
         *
         * @param integers Tags for the integer parameters
         * @param floats Tags for the floating parameters
         */
         ConfigurationFile(const std::string type, const std::vector<std::string>& integers, const std::vector<std::string>& floats);

         /**
         * @brief Destructor
         */
         virtual ~ConfigurationFile() {};

         /**
          * @brief Get integer parameters values
          */
         const ArrayI& integers() const;

         /**
          * @brief Get floating parameters values
          */
         const Array& floats() const;

         /**
          * @brief Read content of configuration file
          */
         virtual void read();

         /**
          * @brief Output run information
          */
         void printInfo() const;
         
      protected:
         /**
          * @brief Resize arrays depending on file type
          */
         void setupStorage();

         /**
          * @brief Read integer values
          */
         void readIntegers();

         /**
          * @brief Read float values
          */
         void readFloats();

      private:
         /**
          * @brief Tags of the integer values
          */
         std::vector<std::string>   mIntegerIDs;

         /**
          * @brief Tags of the float values
          */
         std::vector<std::string>   mFloatIDs;

         /**
          * @brief Storage for the integer values
          */
         ArrayI   mIntegers;

         /**
          * @brief Storage for the float values
          */
         Array   mFloats;

         /**
          * @brief Spread parameters over parallel simulation
          */
         void spreadParameters();
   };

   inline const ArrayI& ConfigurationFile::integers() const
   {
      return this->mIntegers;
   }

   inline const Array& ConfigurationFile::floats() const
   {
      return this->mFloats;
   }

   /// Typedef for a smart pointer of a Parameters
   typedef EPMSHARED_PTR<ConfigurationFile> SmartConfigurationFile;

}

#endif // CONFIGURATIONFILE_HPP
