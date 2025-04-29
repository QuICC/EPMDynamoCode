/** \file HDF5Reader.cpp
 *  \brief Source of the HDF5 reader implementation
 */

// System includes
//
#include <sstream>

// External includes
//

// Class include
//
#include "IO/HDF5/HDF5Reader.hpp"

// Project includes
//
#include "General/EPMException.hpp"
#include "IO/HDF5/HDF5FileDefs.hpp"

namespace EPMDynamo {

   HDF5Reader::HDF5Reader(std::string name, std::string ext, std::string header, std::string version)
      : HDF5File(name, ext, header, version, H5F_ACC_RDONLY), mBlock(0), mMaxCollIORead(0)
   {
   }

   void HDF5Reader::init()
   {
      // Open the file
      this->open();

      // Check file compatibility
      this->checkCompatibility();
   }

   void HDF5Reader::open()
   {
      // Set Access Property list
      this->setFilePList(); 

      // Open file
      this->mFile = H5Fopen(this->filename().c_str(), this->mFlags, this->mFilePList);

      if(this->mFile < 0)
      {
         throw EPMException("HDF5Reader::open", "File opening failed!");
      }
   }

   void HDF5Reader::checkCompatibility()
   {
      // Open the root of the file
      this->mGroup = H5Gopen(this->mFile, "/", H5P_DEFAULT);

      // Try opening the Header attribute
      hid_t  attr;
      attr = H5Aopen(this->mGroup, this->mHeader.c_str(), H5P_DEFAULT);
      if(attr > 0)
      {
         H5Aclose(attr);
      } else
      // If opening failed throw exception
      {
         throw EPMException("HDF5Reader::checkCompatibility", "Wrong file header!");
      }

      // Open version header and compare versions
      attr = H5Aopen(this->mGroup, HDF5FileDefs::VERSIONTAG.c_str(), H5P_DEFAULT);
      float fileVers;
      H5Aread(attr, H5T_NATIVE_FLOAT, &fileVers);
      H5Aclose(attr);
      float runVers;
      std::istringstream conv(this->mVersion);
      conv >> runVers;

      // Throw exception if wrong version found
      if(runVers != fileVers)
      {
         throw EPMException("HDF5Reader::checkCompatibility", "Wrong file version!");

      }

      H5Gclose(this->mGroup);
   }

   void HDF5Reader::finalise()
   {
      this->close();
   }

   void HDF5Reader::close()
   {
      H5Fclose(this->mFile);

      this->freeFilePList();
   }

   void HDF5Reader::readMatrixVector(const std::string dsname, std::vector<MatrixZ> &storage)
   {
      // Set PList to parallel access
      this->setDatasetPList();

      // Set data type correctly
      HDF5TypesConv  genType;
      hid_t  type = genType.atype<EPMComplex>();

      // Compute size of the memory dataspace
      hsize_t iDims[2];

      // Open dataset in file
      this->mDataset = H5Dopen(this->mGroup, dsname.c_str(), H5P_DEFAULT);

      // Get file dataspace
      hid_t  filespace = H5Dget_space(this->mDataset);

      // memory dataspace 
      hid_t  memspace;
      hsize_t memOffset[2];
      memOffset[0] = 0;
      memOffset[1] = 0;

      // Loop over matrices to read from file
      hsize_t pOffset[2];
      pOffset[1] = 0;

      // First do the collective reads
      //
      for(int i = 0; i < this->mMaxCollIORead; ++i)
      {
         // Set full memory space
         iDims[0] = storage.at(i).cols();
         iDims[1] = storage.at(i).rows();
         memspace = H5Screate_simple(2, iDims, NULL);

         // Select memory hyperslabs
         iDims[1] = this->mBlock;
         H5Sselect_hyperslab(memspace, H5S_SELECT_SET, memOffset, NULL, iDims, NULL);

         // Create non regular hyperslab selection in file
         H5Sselect_none(filespace);
         iDims[0] = 1;
         for(unsigned int j =0; j < this->mFileOffsets.at(i).size(); ++j)
         {
            pOffset[0] = this->mFileOffsets.at(i).at(j);
            H5Sselect_hyperslab(filespace, H5S_SELECT_OR, pOffset, NULL, iDims, NULL);
         }

         // Read file hyperslab into memory hyperslab and update offset
         H5Dread(this->mDataset, type, memspace, filespace, this->mDatasetPList, storage.at(i).data());

         // Reset hyperslab to whole dataset
         H5Sselect_all(filespace);

         // Close memory space
         H5Sclose(memspace);
      }

      // Then do the independent reads
      //
      for(unsigned int i = this->mMaxCollIORead; i < this->mFileOffsets.size(); ++i)
      {
         // Set full memory space
         iDims[0] = storage.at(i).cols();
         iDims[1] = storage.at(i).rows();
         memspace = H5Screate_simple(2, iDims, NULL);

         // Select memory hyperslabs
         iDims[1] = this->mBlock;
         H5Sselect_hyperslab(memspace, H5S_SELECT_SET, memOffset, NULL, iDims, NULL);

         // Create non regular hyperslab selection in file
         H5Sselect_none(filespace);
         iDims[0] = 1;
         for(unsigned int j =0; j < this->mFileOffsets.at(i).size(); ++j)
         {
            pOffset[0] = this->mFileOffsets.at(i).at(j);
            H5Sselect_hyperslab(filespace, H5S_SELECT_OR, pOffset, NULL, iDims, NULL);
         }

         // Read file hyperslab into memory hyperslab and update offset
         H5Dread(this->mDataset, type, memspace, filespace, H5P_DEFAULT, storage.at(i).data());

         // Reset hyperslab to whole dataset
         H5Sselect_all(filespace);

         // Close memory space
         H5Sclose(memspace);

      }

      // Close file dataspace
      H5Sclose(filespace);

      // Close Dataset
      H5Dclose(this->mDataset);

      this->freeDatasetPList();
   }
}
