/** \file SourceFileReaderBase.cpp
 *  \brief Source of the base of the source field readers
 */

// System includes
//

#include <iostream>

// External includes
//

// Class include
//
#include "IO/HDF5/Source/SourceFileReaderBase.hpp"

// Project includes
//
#include "IO/HDF5/Source/SourceFileDefs.hpp"

namespace EPMDynamo {

   SourceFileReaderBase::SourceFileReaderBase(std::string name, SmartTruncation pTrunc)
      : SpectralHDF5Reader(SourceFileDefs::BASENAME + name, SourceFileDefs::EXTENSION, SourceFileDefs::HEADER, SourceFileDefs::VERSION, pTrunc)
   {
   }

   void SourceFileReaderBase::readCodensity(std::vector<SpectralSHLShell> &scalar)
   {
      // Open the codensity scalar group
      this->mGroup = H5Gopen(this->mFile, SourceFileDefs::CODENSITYTAG.c_str(), H5P_DEFAULT);

      // Read the codensity expansion
      this->readMatrixVector(SourceFileDefs::CODENSITYTAG, scalar);
      
      // close group
      H5Gclose(this->mGroup);
   }

   void SourceFileReaderBase::readMagnetic(std::vector<SpectralSHLShell> &tor, std::vector<SpectralSHLShell> &pol)
   {
      // Open the magnetic field group
      this->mGroup = H5Gopen(this->mFile, SourceFileDefs::MAGNETICTAG.c_str(), H5P_DEFAULT);

      // Read the Toroidal Magnetic expansion
      this->readMatrixVector(SourceFileDefs::MAGNETICTAG+SourceFileDefs::TOROIDALTAG, tor);

      // Read the Poloidal Magnetic expansion
      this->readMatrixVector(SourceFileDefs::MAGNETICTAG+SourceFileDefs::POLOIDALTAG, pol);
      
      // close group
      H5Gclose(this->mGroup);
   }

   void SourceFileReaderBase::readVelocity(std::vector<SpectralSHLShell> &tor, std::vector<SpectralSHLShell> &pol)
   {
      // Open the velocity field group
      this->mGroup = H5Gopen(this->mFile, SourceFileDefs::VELOCITYTAG.c_str(), H5P_DEFAULT);

      // Read the Toroidal Velocity expansion
      this->readMatrixVector(SourceFileDefs::VELOCITYTAG + SourceFileDefs::TOROIDALTAG, tor);

      // Read the Poloidal Velocity expansion
      this->readMatrixVector(SourceFileDefs::VELOCITYTAG + SourceFileDefs::POLOIDALTAG, pol);
      
      // close group
      H5Gclose(this->mGroup);
   }
}
