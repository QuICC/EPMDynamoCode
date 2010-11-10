/** \file IOBaseSystem.hpp
 *  \brief Implementation of the base properties of the IO control
 */

#ifndef IOBASESYSTEM_HPP
#define IOBASESYSTEM_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "IO/ASCII/StdMessage.hpp"
#include "IO/ASCII/ASCIIWriter.hpp"
#include "IO/ASCII/ParametersFile.hpp"
#include "IO/HDF5/HDF5Writer.hpp"
#include "IO/HDF5/State/StateFileReaderBase.hpp"
#include "IO/HDF5/Imposed/ImposedFieldReaderBase.hpp"
#include "IO/HDF5/Source/SourceFileReaderBase.hpp"
#include "Timestepping/TimestepParameters.hpp"

namespace EPMDynamo {

   /**
    * \brief Implementation of the base properties of the IO control
    */
   class IOBaseSystem
   {
      public:
         /**
         * @brief Constructor
         */
         IOBaseSystem();

         /**
         * @brief Destructor
         */
         virtual ~IOBaseSystem() {};

         /**
          * @brief Get the configuration file data
          */
         const SmartParametersFile cfg() const;

         /**
          * @brief Initialise the base system
          */
         void finaliseBase();

      protected:

      private:
         /**
          * @brief Initialise the base system
          */
         void initBase();

         /**
          * @brief Handle to the configuration file
          */
         SmartParametersFile    mpCfgFile;

         /**
          * @brief Handle to StdMessage buffer
          */
         SmartStdMessage   mpStdMessage;
   };

   inline const SmartParametersFile IOBaseSystem::cfg() const
   {
      return this->mpCfgFile;
   }

}

#endif // IOBASESYSTEM_HPP
