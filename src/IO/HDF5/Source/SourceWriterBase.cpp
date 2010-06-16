/** \file SourceWriterBase.cpp
 *  \brief Source of the base of the source field writer
 */

// System includes
//

#include <iostream>

// External includes
//

// Class include
//
#include "IO/HDF5/SourceWriterBase.hpp"

// Project includes
//
#include "IO/HDF5/SourceDefs.hpp"

namespace EPMDynamo {

   SourceWriterBase::SourceWriterBase(std::string name, SmartSTrunc pSTrunc)
      : SpectralHDF5NWriter(SourceDefs::BASENAME + name + "_NEW", SourceDefs::EXTENSION, SourceDefs::HEADER, SourceDefs::VERSION, pSTrunc)
   {
   }

   void SourceWriterBase::writeCodensity(const std::vector<SpectralSHLShell> &scalar)
   {
      // Create the Codensity scalar group
      this->mGroup = H5Gcreate(this->mFile, SourceDefs::CODENSITYTAG.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

      // Write the codensity expansion
      this->writeMatrixVector(SourceDefs::CODENSITYTAG, scalar);
      
      // close group
      H5Gclose(this->mGroup);
   }

   void SourceWriterBase::writeMagnetic(const std::vector<SpectralSHLShell> &tor, const std::vector<SpectralSHLShell> &pol)
   {
      // Create the Magnetic Field group
      this->mGroup = H5Gcreate(this->mFile, SourceDefs::MAGNETICTAG.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

      // Write the Toroidal Magnetic expansion
      this->writeMatrixVector(SourceDefs::MAGNETICTAG+SourceDefs::TOROIDALTAG, tor);

      // Write the Poloidal Magnetic expansion
      this->writeMatrixVector(SourceDefs::MAGNETICTAG+SourceDefs::POLOIDALTAG, pol);
      
      // close group
      H5Gclose(this->mGroup);
   }

   void SourceWriterBase::writeVelocity(const std::vector<SpectralSHLShell> &tor, const std::vector<SpectralSHLShell> &pol)
   {
      // Create the Velocity Field group
      this->mGroup = H5Gcreate(this->mFile, SourceDefs::VELOCITYTAG.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

      // Write the Toroidal Velocity expansion
      this->writeMatrixVector(SourceDefs::VELOCITYTAG + SourceDefs::TOROIDALTAG, tor);

      // Write the Poloidal Velocity expansion
      this->writeMatrixVector(SourceDefs::VELOCITYTAG + SourceDefs::POLOIDALTAG, pol);
      
      // close group
      H5Gclose(this->mGroup);
   }
}
