/** \file SpectralHDF5NWriter.cpp 
 *  \brief Source of the implementation of a general spectral data file writer
 */

// System includes
//

#include <iostream>

// External includes
//

// Class include
//
#include "IO/HDF5/SpectralHDF5NWriter.hpp"

// Project includes
//
#include "IO/HDF5/SpectralFileDefs.hpp"

namespace EPMDynamo {

   SpectralHDF5NWriter::SpectralHDF5NWriter(std::string name, std::string ext, std::string header, std::string version, SmartTruncation pTrunc)
      : HDF5NWriter(name, ext, header, version), mpTrunc(pTrunc)
   {
      // Set dataset dimensions
      this->setDatasetSize();

      // Set dataset offsets
      this->setDatasetOffsets();

      // This set maximum IO operations
      this->setMaxIO();
   }

   void SpectralHDF5NWriter::setMaxIO()
   {
      this->mMaxCollIOWrite = 0;

      // Get the minimum number of harmonic degrees over all CPUs
      for(int i = 0; i < this->runTrunc()->para().nCore(); ++i)
      {
         if(this->mMaxCollIOWrite < this->runTrunc()->remote(i)->spec()->nL())
         {
            this->mMaxCollIOWrite = this->runTrunc()->remote(i)->spec()->nL();
         }
      }
   }

   void SpectralHDF5NWriter::setDatasetSize()
   {
      // Set the first dimension size
      this->mFileDims.push_back(this->runTrunc()->sim()->hoz()->nH());

      // Set the second dimension size (fastest in C ordering)
      this->mFileDims.push_back(this->runTrunc()->sim()->rad()->nN());
   }

   void SpectralHDF5NWriter::setDatasetOffsets()
   {
      // Get array of local L stored
      ArrayI   ls = this->runTrunc()->local()->spec()->lArray();

      hsize_t  offset;

      // Loop over the stored Ls
      std::vector<hsize_t>  offM;
      for(int i=0; i < ls.size(); ++i)
      {
         // Compute offset from previous harmonic degrees
         offset = 0;
         for(int l = 0; l < ls(i); ++l)
         {
            offset += this->runTrunc()->sim()->hoz()->nM(l);
         }

         // Compute the offset for the stored harmonic orders
         for(int j = 0; j < this->runTrunc()->local()->spec()->nM(i); ++j)
         {
            offM.push_back(offset + this->runTrunc()->local()->spec()->m(j,i));
         }

         // Add offset to vector
         this->mFileOffsets.push_back(offM);
         offM.clear();
      }
   }

   void SpectralHDF5NWriter::writeTruncation(const int N, const int L, const int M, const int Mp)
   {
      // Create the Truncation parameters group
      this->mGroup = H5Gcreate(this->mFile, SpectralFileDefs::TRUNCATIONTAG.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

      // Write Radial truncation value to file
      this->writeScalar(SpectralFileDefs::TRUNCNTAG, N);

      // Write latitudinal truncation value to file
      this->writeScalar(SpectralFileDefs::TRUNCLTAG, L);

      // Write longitudinal truncation value to file
      this->writeScalar(SpectralFileDefs::TRUNCMTAG, M);

      // Write assumed longitudinal symmetry to file
      this->writeScalar(SpectralFileDefs::TRUNCMPTAG, Mp);
      
      // close group
      H5Gclose(this->mGroup);
   }

}
