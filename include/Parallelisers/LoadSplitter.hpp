/** \file LoadSplitter.hpp
 *  \brief Implementation of the work load splitter over the available CPUs
 */

#ifndef LOADSPLITTER_HPP
#define LOADSPLITTER_HPP

// Configuration includes
//
#include "Config/Parallelisation.h"

// System includes
//
#include <vector>

// External includes
//

// Project includes
//
#include "General/EPMTypedefs.hpp"
#include "General/EPMException.hpp"
#include "Domain/SimulationTruncation.hpp"
#include "Domain/CoreTruncation.hpp"
#include "Domain/LocalTruncation.hpp"
#include "Parallelisers/LoadSplitterBase.hpp"

namespace EPMDynamo {

   /**
    * \brief Implementation of the work load splitter over the available CPUs
    *
    * \tparam TSplitAlgo Type of Splitting algorithms
    */
   template <typename TSplitAlgo> class LoadSplitter: public TSplitAlgo
   {
      public:
         /**
          * @brief Constructor
          *
          * @param pSimTrunc Simulation wide truncation information
          * @param nCore Number of cores used
          * @param id ID of the CPU/Core
          */
         LoadSplitter(SmartSimTrunc pSimTrunc, const int nCore, const int id);

         /**
          * @brief Destructor
          */
         virtual ~LoadSplitter() {};

         /**
          * @brief Construct truncation information for all cores
          *
          * @param rpLocal Pointer for local CPU/Core truncation
          * @param rRemote Vector of truncation pointer for remote CPUs/Cores
          */
         void  createTruncations(SmartLocalTrunc &rpLocal, std::vector<SmartCoreTrunc>  &rRemote) const;
         
      protected:

      private:
         /**
          * @brief Check that the splitting was correct
          *
          * @param pLocal Local truncation information
          * @param pRemote Remote truncation information
          */
         void checkSplitting(SmartLocalTrunc pLocal, const std::vector<SmartCoreTrunc>  &pRemote) const;

         /**
          * @brief Check the splittings integrity
          *
          * The goal of this method is to make sure no element gets lost
          *
          * @param pRemote Remote truncation information
          */
         void checkConsistency(const std::vector<SmartCoreTrunc>  &pRemote) const;

         /**
          * @brief Check local CPU load
          *
          * @param pLocal Local truncation information
          */
         void checkLoad(SmartLocalTrunc pLocal) const;
   };

   template <typename TSplitAlgo> LoadSplitter<TSplitAlgo>::LoadSplitter(SmartSimTrunc pSimTrunc, const int nCore, const int id)
      : TSplitAlgo(pSimTrunc, nCore, id)
   {
   }
   
   template <typename TSplitAlgo> void LoadSplitter<TSplitAlgo>::createTruncations(SmartLocalTrunc &rpLocal, std::vector<SmartCoreTrunc>  &rRemote) const
   {
      // RTP data storage
      int rR0;
      int rNr;
      ArrayI rTh0;
      ArrayI rNth;

      // FDSH data storage
      ArrayI fR0;
      ArrayI fNr;
      ArrayI fM;

      // Spec data storage
      ArrayI sL;
      std::vector<ArrayI> sMs;

      // General purpose core truncation pointer
      SmartCoreTrunc    pCore;

      for(int i = 0; i < this->nCore(); ++i)
      {
         // Compute RTP data box splitting
         this->splitRTP(rR0, rNr, rTh0, rNth, i);

         // Compute FDSH data box splitting
         this->splitFDSH(fR0, fNr, fM, i);

         // Compute Spectral data box splitting
         this->splitSpec(sL, sMs, i);

         // Create remote truncation information
         pCore = SmartCoreTrunc(new CoreTruncation(rR0, rNr, rTh0, rNth, fR0, fNr, fM, sL, sMs));

         rRemote.push_back(pCore);

         // Create local truncation information (Might contain some domain specific information
         if(i == this->id())
         {
            rpLocal = SmartLocalTrunc(new LocalTruncation(rR0, rNr, rTh0, rNth, fR0, fNr, fM, sL, sMs));
         }
      }

      // Check the obtained splitting
      this->check(rpLocal, rRemote);
   }

   template <typename TSplitAlgo> void LoadSplitter<TSplitAlgo>::checkSplitting(SmartLocalTrunc pLocal, const std::vector<SmartCoreTrunc>  &pRemote) const
   {
      // Check RTP splitting consistency
      this->checkConsistency(pRemote);

      // Check RTP load splitting
      this->checkLoad(pLocal);
   }

   template <typename TSplitAlgo> void LoadSplitter<TSplitAlgo>::checkConsistency(const std::vector<SmartCoreTrunc>  &pRemote) const
   {
      // Storage for results of consistency checks
      int consistent[3];

      // Check RTP splitting consistency
      consistent[0] = this->checkRTPConsistency(pRemote);

      // Check FDSH splitting consistency
      consistent[1] = this->checkFDSHConsistency(pRemote);

      // Check Spec splitting consistency
      consistent[2] = this->checkSpecConsistency(pRemote);

      // Get the "global" consistency results from MPI code
      #ifdef EPMDYNAMO_MPI
         MPI_Allreduce(MPI_IN_PLACE, consistent, 3, MPI_INT, MPI_PROD, MPI_COMM_WORLD);
      #endif

      // Make it look nice ;)
      if(this->id() == 0)
      {
         std::cout << "*******************************" << std::endl;
         std::cout << "***  Splitting consistency  ***" << std::endl;
         std::cout << "*******************************" << std::endl;
         if(consistent[0])
         {
            std::cout << "***      RTP  : PASSED!     ***" << std::endl;
         } else
         {
            std::cout << "***      RTP  : FAILED!     ***" << std::endl;
         }
         if(consistent[1])
         {
            std::cout << "***      FDSH : PASSED!     ***" << std::endl;
         } else
         {
            std::cout << "***      FDSH : FAILED!     ***" << std::endl;
         }
         if(consistent[2])
         {
            std::cout << "***      SPEC : PASSED!     ***" << std::endl;
         } else
         {
            std::cout << "***      SPEC : FAILED!     ***" << std::endl;
         }
         std::cout << "*******************************" << std::endl;
         std::cout << std::endl;
      }

      if(consistent[0] * consistent[1] * consistent[2] != 1)
      {
         throw EPMException("LoadSplitter::checkConsistency", "Data/Load splitting is not consistent");
      }
   }

   template <typename TSplitAlgo> void LoadSplitter<TSplitAlgo>::checkLoad(SmartLocalTrunc pLocal) const
   {
      // Make it look nice ;)
      if(this->id() == 0)
      {
         std::cout << "**********************************" << std::endl;
         std::cout << "***** Maximal load inbalance *****" << std::endl;
         std::cout << "**********************************" << std::endl;
      }

      // Storage for inbalance results
      double inbalance[3];

      // Check RTP load splitting
      inbalance[0] = this->computeRTPInbalance(pLocal);

      // Check FDSH load splitting
      inbalance[1] = this->computeFDSHInbalance(pLocal);

      // Check Spec load splitting
      inbalance[2] = this->computeSpecInbalance(pLocal);

      // Get the "global" load inbalance for MPI code
      #ifdef EPMDYNAMO_MPI
         MPI_Allreduce(MPI_IN_PLACE, inbalance, 3, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
      #endif

      // Even nicer
      if(this->id() == 0)
      {
         std::cout << "***\t RTP: " << std::ceil(inbalance[0]*100.0) << "%" << std::endl;
         std::cout << "***\t FDSH: " << std::ceil(inbalance[1]*100.0) << "%" << std::endl;
         std::cout << "***\t Spec: " << std::ceil(inbalance[2]*100.0) << "%" << std::endl;
         std::cout << "**********************************" << std::endl;
      }
   }
}

#endif // LOADSPLITTER_HPP
