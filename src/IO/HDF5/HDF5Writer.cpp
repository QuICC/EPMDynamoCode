/** \file HDF5Writer.cpp
 *  \brief Source of the general HDF5 writer
 */

// System includes
//
#include <sstream>

// External includes
//

// Class include
//
#include "IO/HDF5/HDF5Writer.hpp"

// Project includes
//
#include "IO/HDF5/HDF5FileDefs.hpp"

namespace EPMDynamo {

   HDF5Writer::HDF5Writer(std::string name, std::string ext, std::string header, std::string version)
      : HDF5File(name, ext, header, version, H5F_ACC_TRUNC), mMaxCollIOWrite(10000)
   {
   }

   void HDF5Writer::open()
   {
      // Set Create Property list
      this->setFilePList(); 

      // Create file
      this->mFile = H5Fcreate(this->filename().c_str(), this->mFlags, H5P_DEFAULT, this->mFilePList);
   }

   void HDF5Writer::close()
   {
      this->freeFilePList();

      H5Fclose(this->mFile);
   }

   void HDF5Writer::createFileInfo()
   {
      // Select root of the file
      this->mGroup = H5Gopen(this->mFile, "/", H5P_DEFAULT);

      // Define some used variables
      hid_t attr;
      hsize_t dims;
      hid_t  dspace;

      // The header has no data just a name
      dims = 0;
      dspace = H5Screate_simple(0, &dims, NULL);
      attr = H5Acreate(this->mGroup, this->mHeader.c_str(), H5T_NATIVE_INT, dspace, H5P_DEFAULT, H5P_DEFAULT);
      H5Aclose(attr);
      H5Sclose(dspace);

      // The version information is stored as a float
      dims = 1;
      dspace = H5Screate_simple(0, &dims, NULL);

      // Conversion between string and float is done through an istringstream
      std::istringstream  conv(this->mVersion);
      float vers;
      conv >> vers;
      attr = H5Acreate(this->mGroup, HDF5FileDefs::VERSIONTAG.c_str(), H5T_NATIVE_FLOAT, dspace, H5P_DEFAULT, H5P_DEFAULT);
      H5Awrite(attr,H5T_NATIVE_FLOAT, &vers);
      H5Aclose(attr);
      H5Sclose(dspace);

      H5Gclose(this->mGroup);
   }

   void HDF5Writer::writeArray(const std::string dsname, const Array &inArray)
   {
      // Set data type correctly
      HDF5TypesConv  genType;
      hid_t type = genType.ptype<EPMFloat>();

      // Compute size of the dataspace
      hsize_t dims = inArray.size();

      // Create file dataspace
      hid_t filespace = H5Screate_simple(1, &dims, NULL);

      // Create dataset in file
      mDataset = H5Dcreate(this->mGroup, dsname.c_str(), type, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

      // Create memory dataspace 
      hid_t  memspace;

      // Set memory space
      memspace = H5Screate_simple(1, &dims, NULL);

      // Write memory
      if(doesIO())
      {
         H5Dwrite(this->mDataset, type, memspace, filespace, H5P_DEFAULT, inArray.data());
      }
      
      // Close memspace
      H5Sclose(memspace);

      // Close filespace
      H5Sclose(filespace);

      // Close Dataset
      H5Dclose(this->mDataset);
   }

   void HDF5Writer::writeMatrixVector3D(const std::string dsname, const std::vector<Matrix> &storage)
   {
      // Set PList to parallel access
      this->setDatasetPList();

      // Set data type correctly
      HDF5TypesConv  genType;
      hid_t type = genType.ptype<EPMFloat>();

      // Compute size of the dataspace
      hsize_t iDims[2];
      iDims[1] = this->mFileDims.back();

      // Create file dataspace
      hid_t filespace = H5Screate_simple(this->mFileDims.size(), &this->mFileDims.front(), NULL);

      // Create dataset in file
      mDataset = H5Dcreate(this->mGroup, dsname.c_str(), type, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

      // Memory dataspace 
      hid_t  memspace;
      hsize_t pOffsets[3];
      pOffsets[2] = 0;

      // Compute size of the dataspace
      hsize_t dims[3];
      dims[0] = 1;
      dims[2] = this->mFileDims.back();

      // Loop over matrices to store in file
      for(unsigned int i = 0; i < storage.size() ; ++i)
      {
         // Set memory space
         iDims[0] = storage.at(i).cols();
         memspace = H5Screate_simple(2, iDims, NULL);

         pOffsets[0] = this->mFileOffsets.at(i).at(0);
         pOffsets[1] = this->mFileOffsets.at(i).at(1);

         // Select corresponding hyperslab
         dims[1] = storage.at(i).cols();
         H5Sselect_hyperslab(filespace, H5S_SELECT_SET, pOffsets, NULL, dims, NULL);

         // Write memory into hyperslab
         H5Dwrite(this->mDataset, type, memspace, filespace, mDatasetPList, storage.at(i).data());

         // Reset hyperslab to whole dataset
         H5Sselect_all(filespace);
         H5Sclose(memspace);
      }

      // Add some ZERO IO calls to allow for collective writes
      // Set zero memory space
      iDims[0] = 1;
      memspace = H5Screate_simple(1, iDims, NULL);
      H5Sselect_none(memspace);

      // Select zero file space
      H5Sselect_none(filespace);

      for(unsigned int i = 0; i < (this->mMaxCollIOWrite - storage.size()) ; ++i)
      {
         // Write memory into hyperslab
         H5Dwrite(this->mDataset, type, memspace, filespace, mDatasetPList, storage.at(0).data());
      }

      // Close memspace
      H5Sclose(memspace);

      // Close Dataset
      H5Dclose(this->mDataset);

      // Close filespace
      H5Sclose(filespace);

      // Free the dataset PList
      this->freeDatasetPList();
   }

   void HDF5Writer::writeMatrixVector(const std::string dsname, const std::vector<MatrixZ> &storage)
   {
      // Set PList to parallel access
      this->setDatasetPList();

      // Set data type correctly
      HDF5TypesConv  genType;
      hid_t type = genType.atype<EPMComplex>();

      // Compute size of the dataspace
      hsize_t iDims[2];
      iDims[1] = this->mFileDims.back();

      // Create file dataspace
      hid_t filespace = H5Screate_simple(this->mFileDims.size(), &(this->mFileDims.front()), NULL);

      // Create dataset in file
      mDataset = H5Dcreate(this->mGroup, dsname.c_str(), type, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

      // Create memory dataspace 
      hid_t  memspace;
      hsize_t pOffset[2];
      pOffset[1] = 0;

      // Loop over matrices to store in file
      for(unsigned int i = 0; i < storage.size() ; ++i)
      {
         // Set memory space
         iDims[0] = storage.at(i).cols();
         memspace = H5Screate_simple(2, iDims, NULL);

         // Create non regular hyperslab selection in file
         H5Sselect_none(filespace);
         iDims[0] = 1;
         for(unsigned int j =0; j < this->mFileOffsets.at(i).size(); ++j)
         {
            pOffset[0] = this->mFileOffsets.at(i).at(j);
            H5Sselect_hyperslab(filespace, H5S_SELECT_OR, pOffset, NULL, iDims, NULL);
         }

         // Write memory into hyperslab
         H5Dwrite(this->mDataset, type, memspace, filespace, this->mDatasetPList, storage.at(i).data());

         // Reset hyperslab to whole dataset
         H5Sselect_all(filespace);
         H5Sclose(memspace);
      }

      // Add some ZERO IO calls to allow for collective writes
      // Create zero memory space
      memspace = H5Screate_simple(1, iDims, NULL);
      H5Sselect_none(memspace);

      // Create zero file selection space
      H5Sselect_none(filespace);

      // Write zero IO
      for(unsigned int i = 0; i < (this->mMaxCollIOWrite - storage.size()); ++i)
      {
         // Write memory into hyperslab
         H5Dwrite(this->mDataset, type, memspace, filespace, this->mDatasetPList, storage.at(0).data());
      }

      // close memspace
      H5Sclose(memspace);

      // Close Dataset
      H5Dclose(this->mDataset);

      // Close filespace
      H5Sclose(filespace);

      // Free the dataset PList
      this->freeDatasetPList();
   }
}
