/** \file Truncation.cpp
 *  \brief Source of the truncation implementation
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
#include "Parallelisers/SplittingAlgorithms/SerialSplitting.hpp"
#include "Parallelisers/SplittingAlgorithms/TubularSplitting.hpp"
#include "Parallelisers/SplittingAlgorithms/RadialSplitting.hpp"
#include "Parallelisers/SplittingAlgorithms/SHSplitting.hpp"

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
         // Create serial splitter
         LoadSplitter<SerialSplitting>  splitter(this->sim(), 1 ,0);

         // Create truncation
         splitter.createTruncations(this->mpLocal, this->mRemote);
      } else
      {
         try
         {
         #ifdef EPMDYNAMO_SPLIT_SH
            #ifdef EPMDYNAMO_SPLIT_RADIAL
               // Create tubular massive splitter
               LoadSplitter<TubularSplitting>  splitter(this->sim(), this->para().nCore(), this->para().id());

               splitter.createTruncations(this->mpLocal, this->mRemote);
            #else
               // Create Spherical harmonics splitter
               LoadSplitter<SHSplitting>  splitter(this->sim(), this->para().nCore(), this->para().id());

               splitter.createTruncations(this->mpLocal, this->mRemote);
            #endif // EPMDYNAMO_SPLIT_RADIAL
         #else
            #ifdef EPMDYNAMO_SPLIT_RADIAL
               // Create radial splitter
               LoadSplitter<RadialSplitting>  splitter(this->sim(), this->para().nCore(), this->para().id());

               splitter.createTruncations(this->mpLocal, this->mRemote);
            #else
               // This should never happen!
               assert(false);
            #endif // EPMDYNAMO_SPLIT_RADIAL
         #endif // EPMDYNAMO_SPLIT_SH
         }
         catch(EPMException &e)
         {
            e.printStdMessage();
         }
      }
   }
}
