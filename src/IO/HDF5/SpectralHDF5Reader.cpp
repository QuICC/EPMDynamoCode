/** \file SpectralHDF5Reader.cpp
 *  \brief Source of the implementation of a general spectral data file reader
 */

// System includes
//

#include <iostream>

// External includes
//

// Class include
//
#include "IO/HDF5/SpectralHDF5Reader.hpp"

// Project includes
//
#include "General/EPMException.hpp"
#include "IO/HDF5/SpectralFileDefs.hpp"

namespace EPMDynamo {

   SpectralHDF5Reader::SpectralHDF5Reader(std::string name, std::string ext, std::string header, std::string version, SmartTruncation pTrunc)
      : HDF5Reader(name, ext, header, version), mpTrunc(pTrunc)
   {
      // This set maximum IO operations
      this->setMaxIO();
   }

   void SpectralHDF5Reader::setMaxIO()
   {
      this->mMaxCollIORead = 0;
   }

   void SpectralHDF5Reader::setReadArguments()
   {
      // Set block size (i.e. number of contiguous elements)
      this->mBlock = std::min(this->runTrunc()->sim()->rad()->nN(), this->fileTrunc()->rad()->nN());

      // Get array of local L stored
      ArrayI   ls = this->runTrunc()->local()->spec()->lArray();
      int   mp = this->runTrunc()->sim()->hoz()->mp();

      // Get max values from file truncation
      int fileL = this->fileTrunc()->hoz()->maxL();
      int fileM = this->fileTrunc()->hoz()->maxM();

      // Loop over the stored Ls
      int l_;
      int m_;
      std::vector<hsize_t>  offM;
      hsize_t offset;
      for(int i=0; i < ls.size(); ++i)
      {
         l_ = ls(i);
         if(l_ <= fileL)
         {
            // Compute offset from previous harmonic degrees
            offset = 0;
            for(int l = 0; l < l_; ++l)
            {
               offset += this->fileTrunc()->hoz()->nM(l);
            }

            // Compute the offset for the stored harmonic orders
            for(int j = 0; j < this->runTrunc()->local()->spec()->nM(i); ++j)
            {
               m_ = this->runTrunc()->local()->spec()->m(j,i);
               if(mp*m_ <= fileM)
               {
                  offM.push_back(offset + m_);
               }
            }

            // Add offset to vector
            this->mFileOffsets.push_back(offM);
            offM.clear();
         }
      }

      // Get the "global" local minimum for MPI code
      this->mMaxCollIORead = this->mFileOffsets.size();
      #ifdef EPMDYNAMO_MPI
         MPI_Allreduce(MPI_IN_PLACE, &this->mMaxCollIORead, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
      #endif // EPMDYNAMO_MPI
   }

   void SpectralHDF5Reader::readTruncation()
   {
      // Open the truncation parameters group
      this->mGroup = H5Gopen(this->mFile, SpectralFileDefs::TRUNCATIONTAG.c_str(), H5P_DEFAULT);
   
      int maxN;
      int maxL;
      int maxM;
      int mp;
      // Read Radial truncation value to file
      this->readScalar(SpectralFileDefs::TRUNCNTAG, maxN);

      // Read latitudinal truncation value to file
      this->readScalar(SpectralFileDefs::TRUNCLTAG, maxL);

      // Read longitudinal truncation value to file
      this->readScalar(SpectralFileDefs::TRUNCMTAG, maxM);

      // Read assumed longitudinal symmetry value to file
      this->readScalar(SpectralFileDefs::TRUNCMPTAG, mp);
      
      // close group
      H5Gclose(this->mGroup);

      this->mpFileSim = SmartSimTrunc(new SimulationTruncation(maxN, maxN, maxL, maxM, mp));
   }

   void SpectralHDF5Reader::checkTruncation()
   {
      if(this->runTrunc()->sim()->rad()->maxN() != this->fileTrunc()->rad()->maxN())
      {
         std::cout << "Radial truncation doesn't fit" << std::endl;
         if(this->runTrunc()->sim()->rad()->maxN() > this->fileTrunc()->rad()->maxN())
         {
            std::cout << " ---> Zeros have been added!" << std::endl;
         } else
         {
            std::cout << " ---> File data has been truncated!" << std::endl;
         }
         std::cout << std::endl;
      }
      if(this->runTrunc()->sim()->hoz()->maxL() != this->fileTrunc()->hoz()->maxL())
      {
         std::cout << "Latitudinal truncation doesn't fit" << std::endl;
         if(this->runTrunc()->sim()->hoz()->maxL() > this->fileTrunc()->hoz()->maxL())
         {
            std::cout << " ---> Zeros have been added!" << std::endl;
         } else
         {
            std::cout << " ---> File data has been truncated!" << std::endl;
         }
         std::cout << std::endl;
      }
      if(this->runTrunc()->sim()->hoz()->maxM() != this->fileTrunc()->hoz()->maxM())
      {
         std::cout << "Longitudinal truncation doesn't fit" << std::endl;
         if(this->runTrunc()->sim()->hoz()->maxM() > this->fileTrunc()->hoz()->maxM())
         {
            std::cout << " ---> Zeros have been added!" << std::endl;
         } else
         {
            std::cout << " ---> File data has been truncated!" << std::endl;
         }
         std::cout << std::endl;
      }
      if(this->runTrunc()->sim()->hoz()->mp() != this->fileTrunc()->hoz()->mp())
      {
         throw EPMException("SpectralHDF5Reader::checkTruncation", "Incompatible symmetries");
      }
      std::cout << std::endl;
   }
}
