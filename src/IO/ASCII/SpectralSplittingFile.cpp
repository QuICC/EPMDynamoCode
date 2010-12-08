/** \file SpectralSplittingFile.cpp
 *  \brief Source of the spectral space splitting information file
 */

// System includes
//

// External includes
//

// Class include
//
#include "IO/ASCII/SpectralSplittingFile.hpp"

// Project includes
//

namespace EPMDynamo {

   SpectralSplittingFile::SpectralSplittingFile(const SmartTruncation pTrunc)
      : SplittingFile(pTrunc, "Spectral")
   {
   }

   void SpectralSplittingFile::write()
   {
      if(this->doesIO())
      {
         // Do pre write operations
         this->preWrite();

         // Write header information
         this->mFile << "# N=" << this->mpTrunc->sim()->rad()->maxN() << std::endl;
         this->mFile << "# L=" << this->mpTrunc->sim()->hoz()->maxL() << std::endl;
         this->mFile << "# M=" << this->mpTrunc->sim()->hoz()->maxM() << std::endl;
         this->mFile << "L" << "\t\t" << "M" << "\t\t" << "cpuID" << std::endl;

         // Loop over all CPUs
         for(int id=0; id < this->mpTrunc->para().nCore();++id)
         {
            for(int iL=0; iL < this->mpTrunc->remote(id)->spec()->nL();++iL)
            {
               for(int iM=0;iM < this->mpTrunc->remote(id)->spec()->nM(iL); ++iM)
               {
                     this->mFile << this->mpTrunc->remote(id)->spec()->lArray()(iL) << "\t\t";
                     this->mFile << this->mpTrunc->remote(id)->spec()->mArray(iL)(iM) << "\t\t" << id << std::endl;
               }
            }
         }

         // Do post write operations
         this->postWrite();
      }
   }

}
