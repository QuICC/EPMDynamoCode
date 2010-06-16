/** \file SourceReaderBase.cpp
 *  \brief Source of the base of the source field readers
 */

// System includes
//

#include <iostream>

// External includes
//

// Class include
//
#include "IO/HDF5/SourceReaderBase.hpp"

// Project includes
//
#include "IO/HDF5/SourceDefs.hpp"

namespace EPMDynamo {

   SourceReaderBase::SourceReaderBase(std::string name, SmartSTrunc pSTrunc)
      : SpectralHDF5Reader(SourceDefs::BASENAME + name, SourceDefs::EXTENSION, SourceDefs::HEADER, SourceDefs::VERSION, pSTrunc)
   {
   }

   void SourceReaderBase::readCodensity(std::vector<SpectralSHLShell> &scalar)
   {
      // Open the codensity scalar group
      this->mGroup = H5Gopen(this->mFile, SourceDefs::CODENSITYTAG.c_str(), H5P_DEFAULT);

      // Read the codensity expansion
      this->readMatrixVector(SourceDefs::CODENSITYTAG, scalar);
      
      // close group
      H5Gclose(this->mGroup);
   }

   void SourceReaderBase::readMagnetic(std::vector<SpectralSHLShell> &tor, std::vector<SpectralSHLShell> &pol)
   {
      // Open the magnetic field group
      this->mGroup = H5Gopen(this->mFile, SourceDefs::MAGNETICTAG.c_str(), H5P_DEFAULT);

      // Read the Toroidal Magnetic expansion
      this->readMatrixVector(SourceDefs::MAGNETICTAG+SourceDefs::TOROIDALTAG, tor);

      // Read the Poloidal Magnetic expansion
      this->readMatrixVector(SourceDefs::MAGNETICTAG+SourceDefs::POLOIDALTAG, pol);
      
      // close group
      H5Gclose(this->mGroup);
   }

   void SourceReaderBase::readVelocity(std::vector<SpectralSHLShell> &tor, std::vector<SpectralSHLShell> &pol)
   {
      // Open the velocity field group
      this->mGroup = H5Gopen(this->mFile, SourceDefs::VELOCITYTAG.c_str(), H5P_DEFAULT);

      // Read the Toroidal Velocity expansion
      this->readMatrixVector(SourceDefs::VELOCITYTAG + SourceDefs::TOROIDALTAG, tor);

      // Read the Poloidal Velocity expansion
      this->readMatrixVector(SourceDefs::VELOCITYTAG + SourceDefs::POLOIDALTAG, pol);
      
      // close group
      H5Gclose(this->mGroup);
   }
}
