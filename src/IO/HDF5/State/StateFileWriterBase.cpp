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

      if(type == "EQRaRo")
      {
         // Write Ekman scalar to file
         this->writeScalar(StateFileDefs::EKMANTAG, params(0));

         // Write Roberts scalar to file
         this->writeScalar(StateFileDefs::ROBERTSTAG, params(1));

         // Write Rayleigh scalar to file
         this->writeScalar(StateFileDefs::RAYLEIGHTAG, params(2));

         // Write Rossby scalar to file
         this->writeScalar(StateFileDefs::ROSSBYTAG, params(3));
      } else if(type == "EPmPrRa")
      {
         // Write Ekman scalar to file
         this->writeScalar(StateFileDefs::EKMANTAG, params(0));

         // Write magnetic Prandtl scalar to file
         this->writeScalar(StateFileDefs::MAGNETICPRANDTLTAG, params(1));

         // Write Prandtl scalar to file
         this->writeScalar(StateFileDefs::PRANDTLTAG, params(2));

         // Write Rayleigh scalar to file
         this->writeScalar(StateFileDefs::RAYLEIGHTAG, params(3));
      } else if(type == "EPmQRa")
      {
         // Write Ekman scalar to file
         this->writeScalar(StateFileDefs::EKMANTAG, params(0));

         // Write magnetic Prandtl scalar to file
         this->writeScalar(StateFileDefs::MAGNETICPRANDTLTAG, params(1));

         // Write Roberts scalar to file
         this->writeScalar(StateFileDefs::ROBERTSTAG, params(2));

         // Write Rayleigh scalar to file
         this->writeScalar(StateFileDefs::RAYLEIGHTAG, params(3));
      } else if(type == "EEkEmRa")
      {
         // Write Ekman scalar to file
         this->writeScalar(StateFileDefs::EKMANTAG, params(0));

         // Write thermal Ekman scalar to file
         this->writeScalar(StateFileDefs::THERMALEKMANTAG, params(1));

         // Write magnetic Ekman scalar to file
         this->writeScalar(StateFileDefs::MAGNETICEKMANTAG, params(2));

         // Write Rayleigh scalar to file
         this->writeScalar(StateFileDefs::RAYLEIGHTAG, params(3));
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

   void StateFileWriterBase::writeScalarField(const std::string& name, const std::vector<SpectralSHLShell> &scalar)
   {
      // Create the Codensity scalar group
      this->mGroup = H5Gcreate(this->mFile, name.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

      // Write the codensity expansion
      this->writeMatrixVector(name, scalar);
      
      // close group
      H5Gclose(this->mGroup);
   }

   void StateFileWriterBase::writeTorPolField(const std::string& name, const std::vector<SpectralSHLShell> &tor, const std::vector<SpectralSHLShell> &pol)
   {
      // Create the Magnetic Field group
      this->mGroup = H5Gcreate(this->mFile, name.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

      // Write the Toroidal Magnetic expansion
      this->writeMatrixVector(name+StateFileDefs::TOROIDALTAG, tor);

      // Write the Poloidal Magnetic expansion
      this->writeMatrixVector(name+StateFileDefs::POLOIDALTAG, pol);
      
      // close group
      H5Gclose(this->mGroup);
   }
}
