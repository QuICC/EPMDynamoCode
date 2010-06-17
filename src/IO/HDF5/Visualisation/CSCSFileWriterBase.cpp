/** \file CSCSFileWriterBase.cpp
 *  \brief Source of the base of the CSCS visualisation file format writer
 */

// System includes
//

// External includes
//

// Class include
//
#include "IO/HDF5/Visualisation/CSCSFileWriterBase.hpp"

// Project includes
//
#include "IO/HDF5/Visualisation/CSCSFileDefs.hpp"
#include "IO/HDF5/Visualisation/CSCSFileTools.hpp"

namespace EPMDynamo {

   CSCSFileWriterBase::CSCSFileWriterBase(const std::string gridName, SmartTruncation pTrunc)
      : HDF5NWriter(CSCSFileDefs::BASENAME, CSCSFileDefs::EXTENSION, CSCSFileDefs::HEADER, CSCSFileDefs::VERSION), mpTrunc(pTrunc), mGridName(gridName)
   {
      // Set dataset dimensions
      this->setDatasetSize();

      // Set dataset offsets
      this->setDatasetOffsets();

      // This set maximum IO operations
      this->setMaxIO();
   }

   void CSCSFileWriterBase::setMaxIO()
   {
      this->mMaxCollIOWrite = 0;

      // Get the minimum number of radial sheels over all CPUs
      for(int i = 0; i < this->runTrunc()->para().nCore(); ++i)
      {
         if(this->mMaxCollIOWrite < this->runTrunc()->remote(i)->rtp()->nR())
         {
            this->mMaxCollIOWrite = this->runTrunc()->remote(i)->rtp()->nR();
         }
      }
   }

   void CSCSFileWriterBase::setDatasetSize()
   {
      // Set radial dimension
      this->mFileDims.push_back(mpTrunc->sim()->rad()->nR());

      // Set theta dimension inluding the padding
      this->mFileDims.push_back(mpTrunc->sim()->hoz()->nTh());

      // Set phi dimension inluding the padding
      this->mFileDims.push_back(mpTrunc->sim()->hoz()->nPh() + 1);
   }

   void CSCSFileWriterBase::setDatasetOffsets()
   {
      int r0 = this->mpTrunc->local()->rtp()->r0();
      int nR = this->mpTrunc->local()->rtp()->nR();
      int totNth = this->mpTrunc->sim()->hoz()->nTh();
      ArrayI   th0 = this->mpTrunc->local()->rtp()->th0Array();
      ArrayI   nTh = this->mpTrunc->local()->rtp()->nThArray();

      // Loop over the stored Rs
      std::vector<hsize_t>  offR;
      offR.push_back(0);
      offR.push_back(0);
      for(int i=0; i < nR; ++i)
      {
         // Set offset of the radial shell
         offR.at(0) = r0 + i;

         // Set offset of the theta vectors
         if((th0(i)+nTh(i)) == totNth)
         {
            offR.at(1) = 0;
         } else
         {
            offR.at(1) = totNth - (th0(i)+nTh(i)) + 1;
         }

         // add offsets to vector
         this->mFileOffsets.push_back(offR);
      }
   }

   void CSCSFileWriterBase::createBaseGroup()
   {
      // Create the base grid group
      this->mGroup = H5Gcreate(this->mFile, this->mGridName.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

      // close group
      H5Gclose(this->mGroup);
   }

   void CSCSFileWriterBase::writeGrid(const Array& radial, const Array& theta, const Array& phi)
   {
      // Create the Grid group
      std::string groupPath = this->mGridName + "/" + CSCSFileDefs::GRIDTAG;
      this->mGroup = H5Gcreate(this->mFile, groupPath.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

      // Write the radial array to file
      this->writeArray(CSCSFileDefs::RADIALAXISTAG, radial);

      // Write the radial array to file
      this->writeArray(CSCSFileDefs::THETAAXISTAG, theta);

      // Storage for the corrected version
      Array corrected;

      // Write the radial array to file
      CSCSFileTools::padPhiGrid(corrected, phi);
      this->writeArray(CSCSFileDefs::PHIAXISTAG, corrected);
      
      // close group
      H5Gclose(this->mGroup);
   }

   void CSCSFileWriterBase::writeCodensity(const std::vector<SphericalShell> &scalar)
   {
      // Open grid root
      std::string groupPath = this->mGridName + "/";
      this->mGroup = H5Gopen(this->mFile, groupPath.c_str(), H5P_DEFAULT);

      // Make data suitable for file
      std::vector<SphericalShell>   corrected;
      CSCSFileTools::padRTPField(corrected, scalar, mpTrunc->sim()->hoz()->nTh(), mpTrunc->local()->rtp()->th0Array(), mpTrunc->local()->rtp()->nThArray());

      // Write the codensity expansion
      this->writeMatrixVector3D(CSCSFileDefs::CODENSITYTAG, corrected);
      corrected.clear();
      
      // close group
      H5Gclose(this->mGroup);
   }

   void CSCSFileWriterBase::writeMagnetic(const std::vector<SphericalShell> &r, const std::vector<SphericalShell> &theta, const std::vector<SphericalShell> &phi)
   {
      // Create the Magnetic Field group
      std::string groupPath = this->mGridName + "/" + CSCSFileDefs::MAGNETICTAG;
      this->mGroup = H5Gcreate(this->mFile, groupPath.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

      // data storage for corrected field for file
      std::vector<SphericalShell>   corrected;
      CSCSFileTools::padRTPField(corrected, r, mpTrunc->sim()->hoz()->nTh(), mpTrunc->local()->rtp()->th0Array(), mpTrunc->local()->rtp()->nThArray());

      // Write the radial component
      this->writeMatrixVector3D(CSCSFileDefs::MAGNETICTAG+CSCSFileDefs::RADIALTAG, corrected);
      corrected.clear();

      // Make theta data suitable for file
      CSCSFileTools::padRTPField(corrected, theta, mpTrunc->sim()->hoz()->nTh(), mpTrunc->local()->rtp()->th0Array(), mpTrunc->local()->rtp()->nThArray());

      // Write the radial component
      this->writeMatrixVector3D(CSCSFileDefs::MAGNETICTAG+CSCSFileDefs::THETATAG, corrected);
      corrected.clear();

      // Make phi data suitable for file
      CSCSFileTools::padRTPField(corrected, phi, mpTrunc->sim()->hoz()->nTh(), mpTrunc->local()->rtp()->th0Array(), mpTrunc->local()->rtp()->nThArray());

      // Write the radial component
      this->writeMatrixVector3D(CSCSFileDefs::MAGNETICTAG+CSCSFileDefs::PHITAG, corrected);
      corrected.clear();
      
      // close group
      H5Gclose(this->mGroup);
   }

   void CSCSFileWriterBase::writeVelocity(const std::vector<SphericalShell> &r, const std::vector<SphericalShell> &theta, const std::vector<SphericalShell> &phi)
   {
      // Create the Velocity Field group
      std::string groupPath = this->mGridName + "/" + CSCSFileDefs::VELOCITYTAG;
      this->mGroup = H5Gcreate(this->mFile, groupPath.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

      // data storage for corrected field for file
      std::vector<SphericalShell>   corrected;
      CSCSFileTools::padRTPField(corrected, r, mpTrunc->sim()->hoz()->nTh(), mpTrunc->local()->rtp()->th0Array(), mpTrunc->local()->rtp()->nThArray());

      // Write the radial component
      this->writeMatrixVector3D(CSCSFileDefs::VELOCITYTAG+CSCSFileDefs::RADIALTAG, corrected);
      corrected.clear();

      // Write the radial component
      CSCSFileTools::padRTPField(corrected, theta, mpTrunc->sim()->hoz()->nTh(), mpTrunc->local()->rtp()->th0Array(), mpTrunc->local()->rtp()->nThArray());
      this->writeMatrixVector3D(CSCSFileDefs::VELOCITYTAG+CSCSFileDefs::THETATAG, corrected);
      corrected.clear();

      // Write the radial component
      CSCSFileTools::padRTPField(corrected, phi, mpTrunc->sim()->hoz()->nTh(), mpTrunc->local()->rtp()->th0Array(), mpTrunc->local()->rtp()->nThArray());
      this->writeMatrixVector3D(CSCSFileDefs::VELOCITYTAG+CSCSFileDefs::PHITAG, corrected);
      corrected.clear();
      
      // close group
      H5Gclose(this->mGroup);
   }
}
