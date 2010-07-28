/** \file HDF5File.cpp
 *  \brief Source of the implementation of a general HDF5 file
 */

// System includes
//

#include <iostream>

// External includes
//

// Class include
//
#include "IO/HDF5/HDF5File.hpp"

// Project includes
//

namespace EPMDynamo {

   HDF5File::HDF5File(std::string name, std::string ext, std::string header, std::string version, const unsigned int flags)
      : FileBase(name, ext), mHeader(header), mVersion(version), mFlags(flags), mDatasetPList(H5P_DEFAULT), mFilePList(H5P_DEFAULT), mFile(), mGroup(), mDataset()
   {
   }

   void HDF5File::setFilePList()
   {
      #ifdef EPMDYNAMO_MPI
         this->mFilePList = H5Pcreate(H5P_FILE_ACCESS);

         H5Pset_fapl_mpio(this->mFilePList, MPI_COMM_WORLD, MPI_INFO_NULL);
      #endif // EPMDYNAMO_MPI
   }

   void HDF5File::freeFilePList()
   {
      H5Pclose(this->mFilePList);
   }

   void HDF5File::setDatasetPList()
   {
      #ifdef EPMDYNAMO_MPI
         this->mDatasetPList = H5Pcreate(H5P_DATASET_XFER);

         H5Pset_dxpl_mpio(this->mDatasetPList, H5FD_MPIO_COLLECTIVE);
         //H5Pset_dxpl_mpio(this->mDatasetPList, H5FD_MPIO_INDEPENDENT);
      #endif // EPMDYNAMO_MPI
   }

   void HDF5File::freeDatasetPList()
   {
      H5Pclose(this->mDatasetPList);
   }
}
