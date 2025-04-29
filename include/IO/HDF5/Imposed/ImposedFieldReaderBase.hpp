/** \file ImposedFieldReaderBase.hpp
 *  \brief Base of the imposed field readers
 */

#ifndef IMPOSEDFIELDREADERBASE_HPP
#define IMPOSEDFIELDREADERBASE_HPP

// Configuration includes
//
#include "Config/SmartPointer.h"

// System includes
//

// External includes
//

// Project includes
//
#include "General/EPMTypedefs.hpp"
#include "IO/HDF5/SpectralHDF5Reader.hpp"
#include "Domain/Truncation.hpp"

namespace EPMDynamo {

   /**
    * \brief Base of the imposed field readers
    */
   class ImposedFieldReaderBase: public SpectralHDF5Reader
   {
      public:
         /**
         * @brief Constructor
         *
         * @param name File name
         * @param pTrunc Truncation information
         */
         ImposedFieldReaderBase(std::string name, SmartTruncation   pTrunc);

         /**
         * @brief Destructor
         */
         virtual ~ImposedFieldReaderBase() {};

         /**
          * @brief Write State to file
          */
         virtual void read() = 0;
         
      protected:
         /**
          * @brief Read group for Codensity scalar and set mGroup
          *
          * @param scalar Storage to read codensity values
          */
         void readCodensity(std::vector<SpectralSHLShell> &scalar);

         /**
          * @brief Read group for Magnetic field and set mGroup
          *
          * @param tor Storage to read toroidal magnetic values
          * @param pol Storage to read poloidal magnetic values
          */
         void readMagnetic(std::vector<SpectralSHLShell> &tor, std::vector<SpectralSHLShell> &pol);

         /**
          * @brief Read group for Velocity field and set mGroup
          *
          * @param tor Storage to read toroidal velocity values
          * @param pol Storage to read poloidal velocity values
          */
         void readVelocity(std::vector<SpectralSHLShell> &tor, std::vector<SpectralSHLShell> &pol);

      private:
   };

   /// Typedef for a smart pointer of a StdMessage
   typedef EPMSHARED_PTR<ImposedFieldReaderBase> SmartImposedFieldReaderBase;

}

#endif // IMPOSEDFIELDREADERBASE_HPP
