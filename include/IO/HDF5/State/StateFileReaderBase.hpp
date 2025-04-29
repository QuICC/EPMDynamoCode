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
          * @brief Read Physical parameters from file
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
          * @brief Read group for scalar field
          *
          * @param name Name of the scalar field
          * @param scalar Storage for the codensity values
          */
         void readScalarField(const std::string& name, std::vector<SpectralSHLShell> &scalar);

         /**
          * @brief Read group for toroidal/poloidal field
          *
          * @param name Name of the Toroidal/Poloidal field
          * @param tor Storage for the toroidal values
          * @param pol Storage for the poloidal values
          */
         void readTorPolField(const std::string& name, std::vector<SpectralSHLShell> &tor, std::vector<SpectralSHLShell> &pol);

         /**
          * @brief Read group for toroidal component of toroidal/poloidal field
          *
          * @param name Name of the Toroidal/Poloidal field
          * @param tor Storage for the toroidal values
          */
         void readTorField(const std::string& name, std::vector<SpectralSHLShell> &tor);

         /**
          * @brief Read group for poloidal component of toroidal/poloidal field
          *
          * @param name Name of the Toroidal/Poloidal field
          * @param pol Storage for the poloidal values
          */
         void readPolField(const std::string& name, std::vector<SpectralSHLShell> &pol);

         /**
          * @brief Set the given component to zero
          *
          * @param comp Component to zeroify
          */
         void setZero(std::vector<SpectralSHLShell> &comp);

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
