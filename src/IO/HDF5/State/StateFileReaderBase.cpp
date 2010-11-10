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

   void StateFileReaderBase::readPhysical(EPMFloat &E, EPMFloat &q, EPMFloat &Ra, EPMFloat &Ro)
   {
      // Open the physical parameters group
      this->mGroup = H5Gopen(this->mFile, StateFileDefs::PHYSICALTAG.c_str(), H5P_DEFAULT);

      // Read  Ekman scalar from file
      this->readScalar(StateFileDefs::EKMANTAG, E);

      // Read  q scalar from file
      this->readScalar(StateFileDefs::ROBERTSTAG, q);

      // Read  Rayleigh scalar from file
      this->readScalar(StateFileDefs::RAYLEIGHTAG, Ra);

      // Read  Rossby scalar from file
      this->readScalar(StateFileDefs::ROSSBYTAG, Ro);
      
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

   void StateFileReaderBase::readScalarField(const std::string& name, std::vector<SpectralSHLShell> &scalar)
   {
      // Open the codensity scalar group
      this->mGroup = H5Gopen(this->mFile, name.c_str(), H5P_DEFAULT);

      // Read the codensity expansion
      this->readMatrixVector(name, scalar);
      
      // close group
      H5Gclose(this->mGroup);
   }

   void StateFileReaderBase::readTorPolField(const std::string& name, std::vector<SpectralSHLShell> &tor, std::vector<SpectralSHLShell> &pol)
   {
      // Open the magnetic field group
      this->mGroup = H5Gopen(this->mFile, name.c_str(), H5P_DEFAULT);

      // Read the Toroidal Magnetic expansion
      this->readMatrixVector(name+StateFileDefs::TOROIDALTAG, tor);

      // Read the Poloidal Magnetic expansion
      this->readMatrixVector(name+StateFileDefs::POLOIDALTAG, pol);
      
      // close group
      H5Gclose(this->mGroup);
   }

   void StateFileReaderBase::readTorField(const std::string& name, std::vector<SpectralSHLShell> &tor)
   {
      // Open the magnetic field group
      this->mGroup = H5Gopen(this->mFile, name.c_str(), H5P_DEFAULT);

      // Read the Toroidal Magnetic expansion
      this->readMatrixVector(name+StateFileDefs::TOROIDALTAG, tor);
      
      // close group
      H5Gclose(this->mGroup);
   }

   void StateFileReaderBase::readPolField(const std::string& name, std::vector<SpectralSHLShell> &pol)
   {
      // Open the magnetic field group
      this->mGroup = H5Gopen(this->mFile, name.c_str(), H5P_DEFAULT);

      // Read the Poloidal Magnetic expansion
      this->readMatrixVector(name+StateFileDefs::POLOIDALTAG, pol);
      
      // close group
      H5Gclose(this->mGroup);
   }

   void StateFileReaderBase::setZero(std::vector<SpectralSHLShell> &comp)
   {
      // Loop over the shells
      for(unsigned int i=0; i < comp.size(); ++i)
      {
         // Set values to zero
         comp.at(i).setConstant(0.0);
      }
   }
}
