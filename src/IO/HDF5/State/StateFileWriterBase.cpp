/** \file StateFileWriterBase.cpp
 *  \brief Source of the base of the HDF5 state file writer
 */

// System includes
//

#include <iostream>

// External includes
//

// Class include
//
#include "IO/HDF5/State/StateFileWriterBase.hpp"

// Project includes
//
#include "IO/HDF5/State/StateFileDefs.hpp"

namespace EPMDynamo {

   StateFileWriterBase::StateFileWriterBase(SmartTruncation pTrunc)
      : SpectralHDF5NWriter(StateFileDefs::BASENAME, StateFileDefs::EXTENSION, StateFileDefs::HEADER, StateFileDefs::VERSION, pTrunc)
   {
   }

   void StateFileWriterBase::writePhysical(const std::string type, const Array& params)
   {
      // Create the Physical parameters group
      this->mGroup = H5Gcreate(this->mFile, StateFileDefs::PHYSICALTAG.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

      if(type == "EkQRaRo")
      {
         // Write Ekman scalar to file
         this->writeScalar(StateFileDefs::EKMANTAG, params(0));

         // Write q scalar to file
         this->writeScalar(StateFileDefs::ROBERTSTAG, params(1));

         // Write Rayleigh scalar to file
         this->writeScalar(StateFileDefs::RAYLEIGHTAG, params(2));

         // Write Rossby scalar to file
         this->writeScalar(StateFileDefs::ROSSBYTAG, params(3));
      }
      
      // close group
      H5Gclose(this->mGroup);
   }

   void StateFileWriterBase::writeRun(const EPMFloat time, const EPMFloat step)
   {
      // Create the Run parameters group
      this->mGroup = H5Gcreate(this->mFile, StateFileDefs::RUNTAG.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

      // Write reached simulation time to file
      this->writeScalar(StateFileDefs::RUNTIMETAG, time);

      // Write last timestep to file
      this->writeScalar(StateFileDefs::RUNSTEPTAG, step);
      
      // close group
      H5Gclose(this->mGroup);
   }

   void StateFileWriterBase::writeCodensity(const std::vector<SpectralSHLShell> &scalar)
   {
      // Create the Codensity scalar group
      this->mGroup = H5Gcreate(this->mFile, StateFileDefs::CODENSITYTAG.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

      // Write the codensity expansion
      this->writeMatrixVector(StateFileDefs::CODENSITYTAG, scalar);
      
      // close group
      H5Gclose(this->mGroup);
   }

   void StateFileWriterBase::writeMagnetic(const std::vector<SpectralSHLShell> &tor, const std::vector<SpectralSHLShell> &pol)
   {
      // Create the Magnetic Field group
      this->mGroup = H5Gcreate(this->mFile, StateFileDefs::MAGNETICTAG.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

      // Write the Toroidal Magnetic expansion
      this->writeMatrixVector(StateFileDefs::MAGNETICTAG+StateFileDefs::TOROIDALTAG, tor);

      // Write the Poloidal Magnetic expansion
      this->writeMatrixVector(StateFileDefs::MAGNETICTAG+StateFileDefs::POLOIDALTAG, pol);
      
      // close group
      H5Gclose(this->mGroup);
   }

   void StateFileWriterBase::writeVelocity(const std::vector<SpectralSHLShell> &tor, const std::vector<SpectralSHLShell> &pol)
   {
      // Create the Velocity Field group
      this->mGroup = H5Gcreate(this->mFile, StateFileDefs::VELOCITYTAG.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

      // Write the Toroidal Velocity expansion
      this->writeMatrixVector(StateFileDefs::VELOCITYTAG + StateFileDefs::TOROIDALTAG, tor);

      // Write the Poloidal Velocity expansion
      this->writeMatrixVector(StateFileDefs::VELOCITYTAG + StateFileDefs::POLOIDALTAG, pol);
      
      // close group
      H5Gclose(this->mGroup);
   }
}
