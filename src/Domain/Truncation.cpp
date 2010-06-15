/** \file Truncation.cpp
 *  \brief Source of the truncation
 */

// System includes
//

// External includes
//

// Class include
//
#include "Domain/Truncation.hpp"

// Project includes
//
#include "Parallelisers/LoadSplitter.hpp"
#include "Parallelisers/SingleLoadSplitter.hpp"
#include "Parallelisers/FDSHLoadSplitter.hpp"
#include "Parallelisers/MassiveLoadSplitter.hpp"
#include "Parallelisers/TubularLoadSplitter.hpp"

namespace EPMDynamo {

   Truncation::Truncation(const int maxN, const int nR, const int maxL, const int maxM, const int mp, const int nCore)
      : mPara(nCore), mpSimulation(new SimulationTruncation(maxN, nR, maxL, maxM, mp))
   {
      // Split workload
      this->splitLoad();
   }

   void Truncation::splitLoad()
   {
      if(this->para().nCore() == 1)
      {
         LoadSplitter<SingleLoadSplitter>  splitter(this->sim(), 1 ,0);

         splitter.createTruncations(this->mpLocal, this->mRemote);
      } else
      {
         #ifdef EPMDYNAMO_TWOSTEP
            #ifdef EPMDYNAMO_TUBULAR
               LoadSplitter<TubularLoadSplitter>  splitter(this->sim(), this->para().nCore(), this->para().id());
            #else
               LoadSplitter<MassiveLoadSplitter>  splitter(this->sim(), this->para().nCore(), this->para().id());
            #endif // EPMDYNAMO_TUBULAR
         #else
            LoadSplitter<FDSHLoadSplitter>  splitter(this->sim(), this->para().nCore(), this->para().id());
         #endif // EPMDYNAMO_TWOSTEP

         try
         {
            splitter.createTruncations(this->mpLocal, this->mRemote);
         }
         catch(EPMException &e)
         {
            e.printStdMessage();
         }
      }
   }
}
