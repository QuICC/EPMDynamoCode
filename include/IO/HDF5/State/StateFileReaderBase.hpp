/** \file StateFileReaderBase.hpp
 *  \brief Base of the state file readers
 */

#ifndef STATEFILEREADERBASE_HPP
#define STATEFILEREADERBASE_HPP

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
    * \brief Base of the state file readers
    */
   class StateFileReaderBase: public SpectralHDF5Reader
   {
      public:
         /**
         * @brief Constructor
         *
         * @param number Id number of the file to read
         * @param pTrunc Truncation information
         */
         StateFileReaderBase(std::string number, SmartTruncation   pTrunc);

         /**
         * @brief Destructor
         */
         virtual ~StateFileReaderBase() {};

         /**
          * @brief Write State to file
          */
         virtual void read() = 0;

         /**
          * @brief Get state file time
          */
         EPMFloat time() const;

         /**
          * @brief Get state file timestep
          */
         EPMFloat timestep() const;
         
      protected:
         /**
          * @brief Time read from file
          */
         EPMFloat mTime;

         /**
          * @brief Timestep read from file
          */
         EPMFloat mTimestep;

         /**
          * @brief Read Physical parameters to file
          *
          * @param E Storage to read Ekman number
          * @param q Storage to read Roberts number
          * @param Ra Storage to read Rayleigh number
          * @param Ro Storage to read Rossby number
          */
         void readPhysical(EPMFloat &E, EPMFloat &q, EPMFloat &Ra, EPMFloat &Ro);

         /**
          * @brief Read run information to file
          */
         void readRun();

         /**
          * @brief Read group for Codensity scalar and set mGroup
          *
          * @param scalar Storage for the codensity values
          */
         void readCodensity(std::vector<SpectralSHLShell> &scalar);

         /**
          * @brief Read group for Magnetic field and set mGroup
          *
          * @param tor Storage for the magnetic values: toroidal component
          * @param pol Storage for the magnetic values: poloidal component
          */
         void readMagnetic(std::vector<SpectralSHLShell> &tor, std::vector<SpectralSHLShell> &pol);

         /**
          * @brief Read group for Velocity field and set mGroup
          *
          * @param tor Storage for the velocity values: toroidal component
          * @param pol Storage for the velocity values: poloidal component
          */
         void readVelocity(std::vector<SpectralSHLShell> &tor, std::vector<SpectralSHLShell> &pol);

      private:
   };

   inline EPMFloat StateFileReaderBase::time() const
   {
      return this->mTime;
   }

   inline EPMFloat StateFileReaderBase::timestep() const
   {
      return this->mTimestep;
   }

   /// Typedef for a smart pointer of a StdMessage
   typedef EPMSHARED_PTR<StateFileReaderBase> SmartStateFileReaderBase;

}

#endif // STATEFILEREADERBASE_HPP
