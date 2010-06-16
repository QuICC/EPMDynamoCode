/** \file StateFileReaderBase.cpp
 *  \brief Source of the base of the state file readers
 */

// System includes
//

#include <iostream>

// External includes
//

// Class include
//
#include "IO/HDF5/State/StateFileReaderBase.hpp"

// Project includes
//
#include "IO/HDF5/State/StateFileDefs.hpp"

namespace EPMDynamo {

   StateFileReaderBase::StateFileReaderBase(std::string number, SmartTruncation pTrunc)
      : SpectralHDF5Reader(StateFileDefs::BASENAME + number, StateFileDefs::EXTENSION, StateFileDefs::HEADER, StateFileDefs::VERSION, pTrunc), mTime(-1.0), mTimestep(-1.0)
   {
   }

   void StateFileReaderBase::readPhysical(DynamoFloat &E, DynamoFloat &q, DynamoFloat &Ra, DynamoFloat &Ro)
   {
      // Open the physical parameters group
      this->mGroup = H5Gopen(this->mFile, StateFileDefs::PHYSICALTAG.c_str(), H5P_DEFAULT);

      // Read  Ekman scalar from file
      this->readScalar(StateFileDefs::EKMANTAG, E);

      // Read  q scalar from file
      this->readScalar(StateFileDefs::ROBERTSTAG, q);

      // Read  Rayleigh scalar from file
      this->readScalar(StateFileDefs::RAYLEIGHTAG, Ra);

      // Read  Rosby scalar from file
      this->readScalar(StateFileDefs::ROSBYTAG, Ro);
      
      // close group
      H5Gclose(this->mGroup);
   }

   void StateFileReaderBase::readRun()
   {
      // Open the run paramters group
      this->mGroup = H5Gopen(this->mFile, StateFileDefs::RUNTAG.c_str(), H5P_DEFAULT);

      // Read the reached simulation time from file
      this->readScalar(StateFileDefs::RUNTIMETAG, this->mTime);

      // Read the last used timestep from file
      this->readScalar(StateFileDefs::RUNSTEPTAG, this->mTimestep);
      
      // close group
      H5Gclose(this->mGroup);
   }

   void StateFileReaderBase::readCodensity(std::vector<SpectralSHLShell> &scalar)
   {
      // Open the codensity scalar group
      this->mGroup = H5Gopen(this->mFile, StateFileDefs::CODENSITYTAG.c_str(), H5P_DEFAULT);

      // Read the codensity expansion
      this->readMatrixVector(StateFileDefs::CODENSITYTAG, scalar);
      
      // close group
      H5Gclose(this->mGroup);
   }

   void StateFileReaderBase::readMagnetic(std::vector<SpectralSHLShell> &tor, std::vector<SpectralSHLShell> &pol)
   {
      // Open the magnetic field group
      this->mGroup = H5Gopen(this->mFile, StateFileDefs::MAGNETICTAG.c_str(), H5P_DEFAULT);

      // Read the Toroidal Magnetic expansion
      this->readMatrixVector(StateFileDefs::MAGNETICTAG+StateFileDefs::TOROIDALTAG, tor);

      // Read the Poloidal Magnetic expansion
      this->readMatrixVector(StateFileDefs::MAGNETICTAG+StateFileDefs::POLOIDALTAG, pol);
      
      // close group
      H5Gclose(this->mGroup);
   }

   void StateFileReaderBase::readVelocity(std::vector<SpectralSHLShell> &tor, std::vector<SpectralSHLShell> &pol)
   {
      // Open the velocity field group
      this->mGroup = H5Gopen(this->mFile, StateFileDefs::VELOCITYTAG.c_str(), H5P_DEFAULT);

      // Read the Toroidal Velocity expansion
      this->readMatrixVector(StateFileDefs::VELOCITYTAG + StateFileDefs::TOROIDALTAG, tor);

      // Read the Poloidal Velocity expansion
      this->readMatrixVector(StateFileDefs::VELOCITYTAG + StateFileDefs::POLOIDALTAG, pol);
      
      // close group
      H5Gclose(this->mGroup);
   }
}
