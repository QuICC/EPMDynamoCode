/** \file para_DataTransposes.cpp
 *  \brief Parallelisation test for data transposes implemented in the data manipulators
 *
 *  \epmBug Precision test not entirely implemented yet
 */

#define EIGEN_DEFAULT_IO_FORMAT IOFormat(10)

// Configuration includes
//
#include "Config/Parallelisation.h"

// System includes
//
#include <fstream>
#include <iostream>
#include <Eigen/Array>

// TestSuite includes
//

// Project includes
//
#include "General/EPMTypedefs.hpp"
#include "General/EPMException.hpp"
#include "Domain/Truncation.hpp"
#include "Simulations/Types/WSHSimInc.hpp"
#include "Simulations/Types/WSHSimulation.hpp"
#include "Simulations/SimulationConstants.hpp"
#include "Simulations/Traits/SimulationTraits.hpp"

namespace epm = EPMDynamo;

int TRANSFORM_STEPS = 1;

/// Notation simplification typedef for simulation type
typedef epm::WSHSimulation SimulationType;
/// Notation simplification typedef for truncation type
typedef epm::SmartTruncation  SmartTruncation;
/// Notation simplification typedef for spectralSHTransform type
typedef  epm::SimulationTraits<SimulationType>::TransformType SSHTransformType;

void combineRTPTransforms(const int entry, SSHTransformType &sshTrans)
{
   #ifdef EPMDYNAMO_RADIAL_GROUPEDCOMM
      sshTrans.sshManipulator().initiateGroupedBSend(entry);
   #endif // EPMDYNAMO_RADIAL_GROUPEDCOMM

   #ifdef EPMDYNAMO_SH_GROUPEDCOMM
      #ifndef EPMDYNAMO_RADIAL_GROUPEDCOMM
         sshTrans.sshManipulator().initiateGroupedBSend(entry);
      #endif // EPMDYNAMO_RADIAL_GROUPEDCOMM
      sshTrans.shManipulator().initiateGroupedBSend(entry);
   #endif // EPMDYNAMO_SH_GROUPEDCOMM
}

void combineSpectralTransforms(const int entry, SSHTransformType &sshTrans)
{
   #ifdef EPMDYNAMO_SH_GROUPEDCOMM
      #ifndef EPMDYNAMO_RADIAL_GROUPEDCOMM
         sshTrans.sshManipulator().initiateGroupedFSend(entry);
      #endif // EPMDYNAMO_RADIAL_GROUPEDCOMM
      sshTrans.shManipulator().initiateGroupedFSend(entry);
   #endif // EPMDYNAMO_SHGROUPEDCOMM

   #ifdef EPMDYNAMO_RADIAL_GROUPEDCOMM
      sshTrans.sshManipulator().initiateGroupedFSend(entry);
   #endif // EPMDYNAMO_RADIAL_GROUPEDCOMM
}

/**
 * @brief FDSH data transpose test
 */
int runFDSHTransposeTest(epm::SmartTruncation  pTrunc, SSHTransformType &sshTrans)
{
   int status = 0;

   // add Packet size of transform
   sshTrans.sshManipulator().addPacks(3);

   // Get temporary storage
   SSHTransformType::FDSHForwardType &rMTmpQ = sshTrans.sshManipulator().provideFTmp();

   std::ofstream  file;
   if(pTrunc->para().id() == 0)
   {
      file.open("debug_0.txt", std::ios::app);
   }
   if(pTrunc->para().id() == 1)
   {
      file.open("debug_1.txt", std::ios::app);
   }

   for(int m = 0; m < rMTmpQ.nM(); ++m)
   {
      rMTmpQ.rMShell(m).setConstant(m);
      file << "################# m = " << m << " ###################" << std::endl;
      file << rMTmpQ.mshell(m) << std::endl;
   }
   file.close();

   // Transpose the FDSHForwardType values to a suitable order for radial transform
   sshTrans.sshManipulator().send(rMTmpQ);
   sshTrans.sshManipulator().releaseTmp(rMTmpQ);


   // Get temporary storage
   SSHTransformType::FDSHForwardType &rMTmpS = sshTrans.sshManipulator().provideFTmp();
   SSHTransformType::FDSHForwardType &rMTmpT = sshTrans.sshManipulator().provideFTmp();

   for(int m = 0; m < rMTmpQ.nM(); ++m)
   {
      rMTmpS.rMShell(m).setConstant(-42);
      rMTmpT.rMShell(m).setConstant(-42);
   }

   // Transpose the FDSHForwardType values to a suitable order for radial transform
   sshTrans.sshManipulator().send(rMTmpS);
   sshTrans.sshManipulator().releaseTmp(rMTmpS);

   // Transpose the FDSHForwardType values to a suitable order for radial transform
   sshTrans.sshManipulator().send(rMTmpT);
   sshTrans.sshManipulator().releaseTmp(rMTmpT);


   // Get temporary storage and receive data
   SSHTransformType::FDSHBackwardType &rLTmpQ = sshTrans.sshManipulator().provideBTmp();

   sshTrans.sshManipulator().receive(rLTmpQ);

   if(pTrunc->para().id() == 0)
   {
      file.open("debug_0.txt", std::ios::app);
   }
   if(pTrunc->para().id() == 1)
   {
      file.open("debug_1.txt", std::ios::app);
   }

   for(int l = 0; l < rLTmpQ.nL(); ++l)
   {
      file << "################# l = " << l << " ###################" << std::endl;
      file << rLTmpQ.lshell(l) << std::endl;
   }
   file.close();

   // Free temporary storage
   sshTrans.sshManipulator().freeTmp(rLTmpQ);

   // Get temporary storage and receive data
   SSHTransformType::FDSHBackwardType &rLTmpS = sshTrans.sshManipulator().provideBTmp();
   sshTrans.sshManipulator().receive(rLTmpS);

   // Free temporary storage
   sshTrans.sshManipulator().freeTmp(rLTmpS);

   // Get temporary storage and receive data
   SSHTransformType::FDSHBackwardType &rLTmpT = sshTrans.sshManipulator().provideBTmp();
   sshTrans.sshManipulator().receive(rLTmpT);

   // Free temporary storage
   sshTrans.sshManipulator().freeTmp(rLTmpT);

   return status;
}

/**
 * @brief SH data transpose test
 */
int runSHTransposeTest(epm::SmartTruncation  pTrunc, SSHTransformType &sshTrans)
{
   int status = 0;

   // Set the number of packs of communication
   sshTrans.shManipulator().addPacks(1);

   // Get temporary FDSHBackward data
   SSHTransformType::FFTBackwardType  &rBTmp = sshTrans.shManipulator().provideBTmp();

   // Setup values for temporary data
   for(int m =0; m < pTrunc->local()->fdsh()->nM(); ++m)
   {
      epm::ArrayI ms = pTrunc->local()->fdsh()->mArray();
      int n_;
      for(int n =0; n < pTrunc->local()->fdsh()->nR(m); ++n)
      {
         n_ = n+pTrunc->local()->fdsh()->r0(m);
         for(int th =0; th < pTrunc->sim()->hoz()->nTh(); ++th)
         {
            rBTmp.rMShell(m)(th, n).real() = static_cast<epm::EPMFloat>(1000*ms(m)+n_+0.001*th);
            rBTmp.rMShell(m)(th, n).imag() = -static_cast<epm::EPMFloat>(1000*ms(m)+n_+0.001*th);
         }
      }
   }

   // Send data
   sshTrans.shManipulator().send(rBTmp);

   // Get temporary storage for transposed data
   SSHTransformType::FFTForwardType  &rFTmp = sshTrans.shManipulator().provideFTmp();

   // Receive sent data
   sshTrans.shManipulator().receive(rFTmp);

   // prepare transmission in other direction
   sshTrans.shManipulator().addPacks(1);

   // Send data back
   sshTrans.shManipulator().send(rFTmp);

   // Free M order data storage
   sshTrans.shManipulator().releaseTmp(rFTmp);

   // Get second temporary L ordered data
   SSHTransformType::FFTBackwardType  &rBTmp2 = sshTrans.shManipulator().provideBTmp();

   // Receive data sent back
   sshTrans.shManipulator().receive(rBTmp2);

   // Compare sent and received data
   for(int m =0; m < pTrunc->local()->fdsh()->nM(); ++m)
   {
      status += std::ceil((rBTmp2.mshell(m)-rBTmp.mshell(m)).real().maxCoeff());
      status += std::ceil((rBTmp2.mshell(m)-rBTmp.mshell(m)).imag().maxCoeff());
   }

   // Free memory
   sshTrans.shManipulator().releaseTmp(rBTmp);
   sshTrans.shManipulator().freeTmp(rBTmp2);

   return status;
}

void setupTransform(SSHTransformType &sshTrans)
{
   bool state;

   #ifdef EPMDYNAMO_RADIAL_GROUPEDCOMM
      // This is not working yet
      assert(false);

      state = false;
      TRANSFORM_STEPS = 2;
   #else
      state = true;
      TRANSFORM_STEPS = 1;
   #endif // EPMDYNAMO_RADIAL_GROUPEDCOMM

   // Set the number of packs
   sshTrans.sshManipulator().setMaxPacks(4, 4);

   // Desactivate automatic send/recv control
   sshTrans.sshManipulator().setDynamicPacks(state);

   // Activate or desactivate the special communcation entries
   sshTrans.sshManipulator().setDesactivator(state);

   // Finish initialisation of manipulator
   sshTrans.sshManipulator().setup();

   #ifdef EPMDYNAMO_SH_GROUPEDCOMM
      // This is not working yet
      assert(false);

      state = false;
      TRANSFORM_STEPS++;
   #else
      state = true;
      TRANSFORM_STEPS = std::max(TRANSFORM_STEPS, 1);
   #endif // EPMDYNAMO_SH_GROUPEDCOMM

   // Set the number of packs
   sshTrans.shManipulator().setMaxPacks(4, 4);

   // Desactivate automatic send/recv control
   sshTrans.shManipulator().setDynamicPacks(state);

   // Activate or desactivate the special communcation entries
   sshTrans.shManipulator().setDesactivator(state);

   // Finish initialisation of manipulator
   sshTrans.shManipulator().setup();


   // Configure the communication nesting
   #ifdef EPMDYNAMO_SH_GROUPEDCOMM
      // Set the SpectralSH manipulator to be aware of nested grouped communication
      sshTrans.sshManipulator().setInterStageIDs(0, TRANSFORM_STEPS - 1);
      #ifdef EPMDYNAMO_RADIAL_GROUPEDCOMM
         // If its a full massive grouped communication introduce index shift
         sshTrans.shManipulator().setEntryShift(1);
      #endif // EPMDYNAMO_RADIAL_GROUPEDCOMM
   #endif // EPMDYNAMO_SH_GROUPEDCOMM
}

/**
 * @brief Does transform test computation
 */
int runParaTest()
{
   // Set some truncation values
   int maxN = 6;
   int maxL = 24;
   int maxM = 24;
   int Mp = 1;
   int nCore = 1;

   #ifdef EPMDYNAMO_MPI
      // For MPI case the number of CPU is set according to how it's run
      MPI_Comm_size(MPI_COMM_WORLD, &nCore);
   #endif //EPMDYNAMO_MPI

   // Initialise Spectral truncation information
   SmartTruncation   pTrunc = SimulationType::createTrunc(maxN, maxL, maxM, Mp, nCore);

   // Create Transform object
   SSHTransformType    sshTrans(pTrunc);

   // Setup the transform object for parallel versions
   setupTransform(sshTrans);

   int status = 0;

   // Make sure CPUs are synchronized before start of test
   epm::EPMDYNAMO_SYNCHRONIZE;

   // Run Radial split transpose test test
   status += runFDSHTransposeTest(pTrunc, sshTrans);

   // Make sure CPUs are synchronized before start of test
   epm::EPMDYNAMO_SYNCHRONIZE;

   // Run SH split transpose test test
   status += runSHTransposeTest(pTrunc, sshTrans);

   // Make sure CPUs are synchronized before gathering results
   epm::EPMDYNAMO_SYNCHRONIZE;

   // Gather total status from MPI run
   #ifdef EPMDYNAMO_MPI
      // For MPI case the number of CPU is set according to how it's run
      MPI_Allreduce(MPI_IN_PLACE, &status, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
   #endif //EPMDYNAMO_MPI

   return status;
}

/**
 * @brief General main, setting up MPI if required
 *
 * The actual program is in runPrecTest to make sure MPI initialisations
 * are called before anything else end finalization after destruction
 */
int main(int argc, char* argv[])
{
   // Initilise everything that can't be done inside a class
   epm::EPMDYNAMO_INITIALISER;

   int code = -1;

   try
   {
      code = runParaTest();
   }
   catch(epm::EPMException &e)
   {
      e.printStdMessage();
   }

   // Create output to be taken up by CTest
   if(code)
   {
      if(epm::EPMDYNAMO_RANK == 0)
      {
         std::cout << "Failed!" << std::endl;
      }
   } else
   {
      if(epm::EPMDYNAMO_RANK == 0)
      {
         std::cout << "Passed!" << std::endl;
      }
   }

   // Finalise everything that can't be done inside a class
   epm::EPMDYNAMO_FINALIZER;

   return code;
}
