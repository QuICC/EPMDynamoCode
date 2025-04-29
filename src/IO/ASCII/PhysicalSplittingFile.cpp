/** \file PhysicalSplittingFile.cpp
 *  \brief Source of the physical space splitting information file
 */

// System includes
//

// External includes
//

// Class include
//
#include "IO/ASCII/PhysicalSplittingFile.hpp"

// Project includes
//

namespace EPMDynamo {

   PhysicalSplittingFile::PhysicalSplittingFile(const SmartTruncation pTrunc)
      : SplittingFile(pTrunc, "Physical")
   {
   }

   void PhysicalSplittingFile::write()
   {
      if(this->doesIO())
      {
         // Do pre write operations
         this->preWrite();

         // Write header information
         this->mFile << "# Nr=" << this->mpTrunc->sim()->rad()->nR() << std::endl;
         this->mFile << "# Nth=" << this->mpTrunc->sim()->hoz()->nTh() << std::endl;
         this->mFile << "# Nph=" << this->mpTrunc->sim()->hoz()->nPh() << std::endl;
         this->mFile << "Nr" << "\t\t" << "Nth" << "\t\t" << "cpuID" << std::endl;

         // Loop over all CPUs
         for(int id=0; id < this->mpTrunc->para().nCore();++id)
         {
            for(int n=0;n < this->mpTrunc->remote(id)->rtp()->nR(); ++n)
            {
               for(int th=0; th <this->mpTrunc->remote(id)->rtp()->nTh(n);++th)
               {
                     this->mFile << this->mpTrunc->remote(id)->rtp()->r0()+n << "\t\t";
                     this->mFile << this->mpTrunc->remote(id)->rtp()->th0(n)+th << "\t\t" << id << std::endl;
               }
            }
         }

         // Do post write operations
         this->postWrite();
      }
   }

}
