/** \file ImposedFieldWriterBase.cpp
 *  \brief Source of the base of the imposed field writers
 */

// System includes
//

#include <iostream>

// External includes
//

// Class include
//
#include "IO/HDF5/Imposed/ImposedFieldWriterBase.hpp"

// Project includes
//
#include "IO/HDF5/Imposed/ImposedFieldDefs.hpp"

namespace EPMDynamo {

   ImposedFieldWriterBase::ImposedFieldWriterBase(std::string name, SmartTruncation pTrunc)
      : SpectralHDF5NWriter(ImposedFieldDefs::BASENAME + name + "_NEW", ImposedFieldDefs::EXTENSION, ImposedFieldDefs::HEADER, ImposedFieldDefs::VERSION, pTrunc)
   {
   }

   void ImposedFieldWriterBase::writeCodensity(const std::vector<SpectralSHLShell> &scalar)
   {
      // Create the Codensity scalar group
      this->mGroup = H5Gcreate(this->mFile, ImposedFieldDefs::CODENSITYTAG.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

      // Write the codensity expansion
      this->writeMatrixVector(ImposedFieldDefs::CODENSITYTAG, scalar);
      
      // close group
      H5Gclose(this->mGroup);
   }

   void ImposedFieldWriterBase::writeMagnetic(const std::vector<SpectralSHLShell> &tor, const std::vector<SpectralSHLShell> &pol)
   {
      // Create the Magnetic Field group
      this->mGroup = H5Gcreate(this->mFile, ImposedFieldDefs::MAGNETICTAG.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

      // Write the Toroidal Magnetic expansion
      this->writeMatrixVector(ImposedFieldDefs::MAGNETICTAG+ImposedFieldDefs::TOROIDALTAG, tor);

      // Write the Poloidal Magnetic expansion
      this->writeMatrixVector(ImposedFieldDefs::MAGNETICTAG+ImposedFieldDefs::POLOIDALTAG, pol);
      
      // close group
      H5Gclose(this->mGroup);
   }

   void ImposedFieldWriterBase::writeVelocity(const std::vector<SpectralSHLShell> &tor, const std::vector<SpectralSHLShell> &pol)
   {
      // Create the Velocity Field group
      this->mGroup = H5Gcreate(this->mFile, ImposedFieldDefs::VELOCITYTAG.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

      // Write the Toroidal Velocity expansion
      this->writeMatrixVector(ImposedFieldDefs::VELOCITYTAG + ImposedFieldDefs::TOROIDALTAG, tor);

      // Write the Poloidal Velocity expansion
      this->writeMatrixVector(ImposedFieldDefs::VELOCITYTAG + ImposedFieldDefs::POLOIDALTAG, pol);
      
      // close group
      H5Gclose(this->mGroup);
   }
}
