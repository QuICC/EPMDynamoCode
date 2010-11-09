/** \file IOSystemBase.hpp
 *  \brief Implementation of the base properties of the IO control
 */

#ifndef IOSYSTEMBASE_HPP
#define IOSYSTEMBASE_HPP

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
   class IOSystemBase
   {
      public:
         /**
         * @brief Constructor
         */
         IOSystemBase();

         /**
         * @brief Destructor
         */
         virtual ~IOSystemBase() {};

         /**
          * @brief Get the configuration file data
          */
         const SmartParametersFile cfg() const;

      protected:

         /**
          * @brief Initialise the base system
          */
         void finaliseBase();

      private:
         /**
          * @brief Initialise the base system
          */
         void init();

         /**
          * @brief Handle to the configuration file
          */
         SmartParametersFile    mpCfgFile;

         /**
          * @brief Handle to StdMessage buffer
          */
         SmartStdMessage   mpStdMessage;
   };

   inline const SmartParametersFile IOSystemBase::cfg() const
   {
      return this->mpCfgFile;
   }

}

#endif // IOSYSTEMBASE_HPP
