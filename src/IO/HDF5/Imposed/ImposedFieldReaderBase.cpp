/** \file ImposedFieldReaderBase.cpp
 *  \brief Source of the base of the imposed field reader
 */

// System includes
//

#include <iostream>

// External includes
//

// Class include
//
#include "IO/HDF5/Imposed/ImposedFieldReaderBase.hpp"

// Project includes
//
#include "IO/HDF5/Imposed/ImposedFieldDefs.hpp"

namespace EPMDynamo {

   ImposedFieldReaderBase::ImposedFieldReaderBase(std::string name, SmartTruncation pTrunc)
      : SpectralHDF5Reader(ImposedFieldDefs::BASENAME + name, ImposedFieldDefs::EXTENSION, ImposedFieldDefs::HEADER, ImposedFieldDefs::VERSION, pTrunc)
   {
   }

   void ImposedFieldReaderBase::readCodensity(std::vector<SpectralSHLShell> &scalar)
   {
      // Open the codensity scalar group
      this->mGroup = H5Gopen(this->mFile, ImposedFieldDefs::CODENSITYTAG.c_str(), H5P_DEFAULT);

      // Read the codensity expansion
      this->readMatrixVector(ImposedFieldDefs::CODENSITYTAG, scalar);
      
      // close group
      H5Gclose(this->mGroup);
   }

   void ImposedFieldReaderBase::readMagnetic(std::vector<SpectralSHLShell> &tor, std::vector<SpectralSHLShell> &pol)
   {
      // Open the magnetic field group
      this->mGroup = H5Gopen(this->mFile, ImposedFieldDefs::MAGNETICTAG.c_str(), H5P_DEFAULT);

      // Read the Toroidal Magnetic expansion
      this->readMatrixVector(ImposedFieldDefs::MAGNETICTAG+ImposedFieldDefs::TOROIDALTAG, tor);

      // Read the Poloidal Magnetic expansion
      this->readMatrixVector(ImposedFieldDefs::MAGNETICTAG+ImposedFieldDefs::POLOIDALTAG, pol);
      
      // close group
      H5Gclose(this->mGroup);
   }

   void ImposedFieldReaderBase::readVelocity(std::vector<SpectralSHLShell> &tor, std::vector<SpectralSHLShell> &pol)
   {
      // Open the velocity field group
      this->mGroup = H5Gopen(this->mFile, ImposedFieldDefs::VELOCITYTAG.c_str(), H5P_DEFAULT);

      // Read the Toroidal Velocity expansion
      this->readMatrixVector(ImposedFieldDefs::VELOCITYTAG + ImposedFieldDefs::TOROIDALTAG, tor);

      // Read the Poloidal Velocity expansion
      this->readMatrixVector(ImposedFieldDefs::VELOCITYTAG + ImposedFieldDefs::POLOIDALTAG, pol);
      
      // close group
      H5Gclose(this->mGroup);
   }
}
