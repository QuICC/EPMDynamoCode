/** \file SourceFileWriterBase.cpp
 *  \brief Source of the base of the source field writer
 */

// System includes
//

#include <iostream>

// External includes
//

// Class include
//
#include "IO/HDF5/Source/SourceFileWriterBase.hpp"

// Project includes
//
#include "IO/HDF5/Source/SourceFileDefs.hpp"

namespace EPMDynamo {

   SourceFileWriterBase::SourceFileWriterBase(std::string name, SmartTruncation pTrunc)
      : SpectralHDF5NWriter(SourceFileDefs::BASENAME + name + "_NEW", SourceFileDefs::EXTENSION, SourceFileDefs::HEADER, SourceFileDefs::VERSION, pTrunc)
   {
   }

   void SourceFileWriterBase::writeCodensity(const std::vector<SpectralSHLShell> &scalar)
   {
      // Create the Codensity scalar group
      this->mGroup = H5Gcreate(this->mFile, SourceFileDefs::CODENSITYTAG.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

      // Write the codensity expansion
      this->writeMatrixVector(SourceFileDefs::CODENSITYTAG, scalar);
      
      // close group
      H5Gclose(this->mGroup);
   }

   void SourceFileWriterBase::writeMagnetic(const std::vector<SpectralSHLShell> &tor, const std::vector<SpectralSHLShell> &pol)
   {
      // Create the Magnetic Field group
      this->mGroup = H5Gcreate(this->mFile, SourceFileDefs::MAGNETICTAG.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

      // Write the Toroidal Magnetic expansion
      this->writeMatrixVector(SourceFileDefs::MAGNETICTAG+SourceFileDefs::TOROIDALTAG, tor);

      // Write the Poloidal Magnetic expansion
      this->writeMatrixVector(SourceFileDefs::MAGNETICTAG+SourceFileDefs::POLOIDALTAG, pol);
      
      // close group
      H5Gclose(this->mGroup);
   }

   void SourceFileWriterBase::writeVelocity(const std::vector<SpectralSHLShell> &tor, const std::vector<SpectralSHLShell> &pol)
   {
      // Create the Velocity Field group
      this->mGroup = H5Gcreate(this->mFile, SourceFileDefs::VELOCITYTAG.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

      // Write the Toroidal Velocity expansion
      this->writeMatrixVector(SourceFileDefs::VELOCITYTAG + SourceFileDefs::TOROIDALTAG, tor);

      // Write the Poloidal Velocity expansion
      this->writeMatrixVector(SourceFileDefs::VELOCITYTAG + SourceFileDefs::POLOIDALTAG, pol);
      
      // close group
      H5Gclose(this->mGroup);
   }
}
