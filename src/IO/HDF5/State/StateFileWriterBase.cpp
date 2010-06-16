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

   void StateFileWriterBase::writePhysical(const DynamoFloat E, const DynamoFloat q, const DynamoFloat Ra, const DynamoFloat Ro)
   {
      // Create the Physical parameters group
      this->mGroup = H5Gcreate(this->mFile, StateFileDefs::PHYSICALTAG.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

      // Write Ekman scalar to file
      this->writeScalar(StateFileDefs::EKMANTAG, E);

      // Write q scalar to file
      this->writeScalar(StateFileDefs::ROBERTSTAG, q);

      // Write Rayleigh scalar to file
      this->writeScalar(StateFileDefs::RAYLEIGHTAG, Ra);

      // Write Rosby scalar to file
      this->writeScalar(StateFileDefs::ROSBYTAG, Ro);
      
      // close group
      H5Gclose(this->mGroup);
   }

   void StateFileWriterBase::writeRun(const DynamoFloat time, const DynamoFloat step)
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
