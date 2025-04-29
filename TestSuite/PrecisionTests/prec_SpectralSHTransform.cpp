/** \file prec_SpectralSHTransform.cpp
 *  \brief Precision test for SpectralSH transforms
 *
 *  \epmTSTodo Precision test not entirely implemented yet
 */

#define EIGEN_DEFAULT_IO_FORMAT IOFormat(10)

// Configuration includes
//
#include "Config/Parallelisation.h"

// System includes
//
#include <iostream>

// TestSuite includes
//
#include "ErrorComputer.hpp"
#include "CartesianExpressions.hpp"

// Project includes
//
#include "General/EPMTypedefs.hpp"
#include "General/EPMException.hpp"
#include "General/MathConstants.hpp"
#include "Domain/Truncation.hpp"
#include "Polynomials/AssocLegendreBasis.hpp"
#include "GeneralScalars/RTPScalar.hpp"
#include "Simulations/Types/WSHSimInc.hpp"
#include "Simulations/Types/WSHSimulation.hpp"
#include "Simulations/SimulationConstants.hpp"
#include "Simulations/Traits/SimulationTraits.hpp"

namespace epm = EPMDynamo;

/// Notation simplification typedef for simulation type
typedef epm::WSHSimulation SimulationType;
/// Notation simplification typedef for truncation type
typedef epm::SmartTruncation  SmartTruncation;
/// Notation simplification typedef for radial basis type
typedef SimulationType::RadialBasisType RadialBasisType;
/// Notation simplification typedef for spectralSHTransform type
typedef  epm::SimulationTraits<SimulationType>::TransformType SSHTransformType;
/// Notation simplification typedef for Spectral scalar type
typedef  epm::WSHSimulation::ScalarType SpectralScalarType;
/// Notation simplification typedef for Spectral tor/pol field type
typedef  epm::TorPolField<epm::WSHSimulation> TorPolFieldType;

/// Counter for the number of performed tests
int PERFORMED_TESTS = 0;

/// Number of steps required to perform the transform
int TRANSFORM_STEPS = 1;

/**
 * @brief initialise RTP scalar
 */
void initRTPScalar(epm::RTPScalar &rRTPValues, SmartTruncation pTrunc)
{
   epm::CartesianExpressions   helper(pTrunc);

   for(int n = 0; n < rRTPValues.nR(); n++)
   {
      rRTPValues.rShell(n).setConstant(0.0);
      rRTPValues.rShell(n) = 3.0*(helper.x(n,2) - helper.y(n,2));
   }
}

/**
 * @brief initialise RTP scalar for gradient computation
 */
void initRTPGradientScalar(epm::RTPScalar &rRTPValues, SmartTruncation pTrunc)
{
   epm::CartesianExpressions   helper(pTrunc);

   for(int n = 0; n < rRTPValues.nR(); n++)
   {
      rRTPValues.rShell(n).setConstant(0.0);
      rRTPValues.rShell(n) = 1.4242*helper.z(n,1) + helper.x(n,3).array()*helper.y(n,1).array();
   }
}

/**
 * @brief initialise RTP Field solution to gradient computation
 */
void initRTPGradientField(epm::RTPField &rRTPField, SmartTruncation pTrunc)
{
   epm::CartesianExpressions   helper(pTrunc);

   int rows = rRTPField.nPh();
   int cols;

   epm::Matrix xComp;
   epm::Matrix yComp;
   epm::Matrix zComp;

   for(int n = 0; n < rRTPField.nR(); n++)
   {
      cols = rRTPField.nTh(n);
      xComp.resize(rows, cols);
      yComp.resize(rows, cols);
      zComp.resize(rows, cols);

      rRTPField.rR().rShell(n).setConstant(0.0);
      rRTPField.rTheta().rShell(n).setConstant(0.0);
      rRTPField.rPhi().rShell(n).setConstant(0.0);

      xComp = 3.0*(helper.x(n, 2).array()*helper.y(n,1)).array();
      yComp = helper.x(n, 3);
      zComp.setConstant(1.4242);

      rRTPField.rR().rShell(n) = helper.xyzToR(xComp, yComp, zComp, n);
      rRTPField.rTheta().rShell(n) = helper.xyzToTh(xComp, yComp, zComp, n);
      rRTPField.rPhi().rShell(n) = helper.xyzToPh(xComp, yComp, zComp, n);
   }
}

/**
 * @brief init SpectralSH scalar 
 */
void initSpectralSHScalar(SpectralScalarType &rSSHValues)
{
   int nN = rSSHValues.trunc()->sim()->rad()->nN();

   epm::ArrayI ls = rSSHValues.trunc()->local()->spec()->lArray();
   for(int l=0; l < rSSHValues.trunc()->local()->spec()->nL(); ++l)
   {
      rSSHValues.rLShell(l).setConstant(epm::EPMComplex(ls(l), -ls(l)));

      if(rSSHValues.trunc()->local()->spec()->mArray(l)(0) == 0)
      {
         for(int n =0; n< nN; ++n)
         {
            rSSHValues.rLShell(l)(n,0).imag() = 0.0;
         }
      }
   }
}

/**
 * @brief initialise RTP Field
 */
void initRTPField(epm::RTPField &rRTPField, SmartTruncation pTrunc)
{
   epm::CartesianExpressions   helper(pTrunc);

   int rows = rRTPField.nPh();
   int cols;

   epm::Matrix xComp;
   epm::Matrix yComp;
   epm::Matrix zComp;

   for(int n = 0; n < rRTPField.nR(); n++)
   {
      cols = rRTPField.nTh(n);
      xComp.resize(rows, cols);
      yComp.resize(rows, cols);
      zComp.resize(rows, cols);

      rRTPField.rR().rShell(n).setConstant(0.0);
      rRTPField.rTheta().rShell(n).setConstant(0.0);
      rRTPField.rPhi().rShell(n).setConstant(0.0);

      xComp = helper.x(n, 2) + helper.y(n,1) + helper.z(n,3);
      yComp = helper.x(n, 2).array()* helper.y(n,1).array();
      zComp = helper.z(n,3);

      rRTPField.rR().rShell(n) = helper.xyzToR(xComp, yComp, zComp, n);
      rRTPField.rTheta().rShell(n) = helper.xyzToTh(xComp, yComp, zComp, n);
      rRTPField.rPhi().rShell(n) = helper.xyzToPh(xComp, yComp, zComp, n);
   }
}

/**
 * @brief initialise Divergence free RTP Field
 */
void initDiv0RTPField(epm::RTPField &rRTPField, SmartTruncation pTrunc)
{
   epm::CartesianExpressions   helper(pTrunc);

   int rows = rRTPField.nPh();
   int cols;

   epm::Matrix xComp;
   epm::Matrix yComp;
   epm::Matrix zComp;

   for(int n = 0; n < rRTPField.nR(); n++)
   {
      cols = rRTPField.nTh(n);
      xComp.resize(rows, cols);
      yComp.resize(rows, cols);
      zComp.resize(rows, cols);

      rRTPField.rR().rShell(n).setConstant(0.0);
      rRTPField.rTheta().rShell(n).setConstant(0.0);
      rRTPField.rPhi().rShell(n).setConstant(0.0);

      xComp = helper.x(n, 1) + helper.y(n,1);
      yComp = -helper.y(n, 1) + helper.z(n,1);
      zComp = helper.x(n,3).array()*helper.y(n,3).array();

      rRTPField.rR().rShell(n) = helper.xyzToR(xComp, yComp, zComp, n);
      rRTPField.rTheta().rShell(n) = helper.xyzToTh(xComp, yComp, zComp, n);
      rRTPField.rPhi().rShell(n) = helper.xyzToPh(xComp, yComp, zComp, n);
   }
}

/**
 * @brief initialise Divergence free RTP Field
 */
void initRTPCurlTorPol(epm::RTPField &rRTPField, SmartTruncation pTrunc)
{
   epm::CartesianExpressions   helper(pTrunc);

   int rows = rRTPField.nPh();
   int cols;

   epm::Matrix xComp;
   epm::Matrix yComp;
   epm::Matrix zComp;

   for(int n = 0; n < rRTPField.nR(); n++)
   {
      cols = rRTPField.nTh(n);
      xComp.resize(rows, cols);
      yComp.resize(rows, cols);
      zComp.resize(rows, cols);

      rRTPField.rR().rShell(n).setConstant(0.0);
      rRTPField.rTheta().rShell(n).setConstant(0.0);
      rRTPField.rPhi().rShell(n).setConstant(0.0);

      xComp = helper.x(n, 1) + helper.y(n,1);
      yComp = -helper.y(n, 1) + helper.z(n,1);
      zComp = helper.x(n,3).array()*helper.y(n,3).array();

      rRTPField.rR().rShell(n) = helper.xyzToR(xComp, yComp, zComp, n);
      rRTPField.rTheta().rShell(n) = helper.xyzToTh(xComp, yComp, zComp, n);
      rRTPField.rPhi().rShell(n) = helper.xyzToPh(xComp, yComp, zComp, n);
   }
}

/**
 * @brief initialise Divergence free RTP Field
 */
void initRTPCurlResult(epm::RTPField &rRTPField, SmartTruncation pTrunc)
{
   epm::CartesianExpressions   helper(pTrunc);

   int rows = rRTPField.nPh();
   int cols;

   epm::Matrix xComp;
   epm::Matrix yComp;
   epm::Matrix zComp;

   for(int n = 0; n < rRTPField.nR(); n++)
   {
      cols = rRTPField.nTh(n);
      xComp.resize(rows, cols);
      yComp.resize(rows, cols);
      zComp.resize(rows, cols);

      rRTPField.rR().rShell(n).setConstant(0.0);
      rRTPField.rTheta().rShell(n).setConstant(0.0);
      rRTPField.rPhi().rShell(n).setConstant(0.0);

      xComp.setConstant(-1);
      xComp += 3.0*(helper.x(n, 3).array()*helper.y(n,2).array());
      yComp = -3.0*(helper.x(n, 2).array()*helper.y(n,3).array());
      zComp.setConstant(-1);

      rRTPField.rR().rShell(n) = helper.xyzToR(xComp, yComp, zComp, n);
      rRTPField.rTheta().rShell(n) = helper.xyzToTh(xComp, yComp, zComp, n);
      rRTPField.rPhi().rShell(n) = helper.xyzToPh(xComp, yComp, zComp, n);
   }
}

/**
 * @brief initialise Divergence free RTP Field
 */
void initRTPNTerms(epm::RTPField &rRTPField, SmartTruncation pTrunc)
{
   epm::CartesianExpressions   helper(pTrunc);

   int rows = rRTPField.nPh();
   int cols;

   epm::Matrix xComp;
   epm::Matrix yComp;
   epm::Matrix zComp;

   for(int n = 0; n < rRTPField.nR(); n++)
   {
      cols = rRTPField.nTh(n);
      xComp.resize(rows, cols);
      yComp.resize(rows, cols);
      zComp.resize(rows, cols);

      rRTPField.rR().rShell(n).setConstant(0.0);
      rRTPField.rTheta().rShell(n).setConstant(0.0);
      rRTPField.rPhi().rShell(n).setConstant(0.0);

      xComp = helper.x(n, 1) + helper.y(n, 3) + helper.z(n, 1) + (helper.x(n, 1).array()* helper.y(n, 1).array()*helper.z(n, 1).array();
      yComp = helper.z(n, 4) + helper.x(n, 1).array() * helper.y(n, 2).array();
      zComp = helper.x(n, 1);

      rRTPField.rR().rShell(n) = helper.xyzToR(xComp, yComp, zComp, n);
      rRTPField.rTheta().rShell(n) = helper.xyzToTh(xComp, yComp, zComp, n);
      rRTPField.rPhi().rShell(n) = helper.xyzToPh(xComp, yComp, zComp, n);
   }
}

/**
 * @brief initialise Divergence free RTP Field
 */
void initRTPNTermsCurlResult(epm::RTPField &rRTPField, SmartTruncation pTrunc)
{
   epm::CartesianExpressions   helper(pTrunc);

   int rows = rRTPField.nPh();
   int cols;

   epm::Matrix xComp;
   epm::Matrix yComp;
   epm::Matrix zComp;

   for(int n = 0; n < rRTPField.nR(); n++)
   {
      cols = rRTPField.nTh(n);
      xComp.resize(rows, cols);
      yComp.resize(rows, cols);
      zComp.resize(rows, cols);

      rRTPField.rR().rShell(n).setConstant(0.0);
      rRTPField.rTheta().rShell(n).setConstant(0.0);
      rRTPField.rPhi().rShell(n).setConstant(0.0);

      xComp = -4.0*helper.z(n, 3);
      yComp = helper.x(n, 1).array() * helper.y(n, 1).array();
      zComp = -2.0*helper.y(n, 2) - helper.x(n, 1).array() * helper.z(n, 1).array();

      rRTPField.rR().rShell(n) = helper.xyzToR(xComp, yComp, zComp, n);
      rRTPField.rTheta().rShell(n) = helper.xyzToTh(xComp, yComp, zComp, n);
      rRTPField.rPhi().rShell(n) = helper.xyzToPh(xComp, yComp, zComp, n);
   }
}

/**
 * @brief initialise Divergence free RTP Field
 */
void initRTPNTermsCurlCurlResult(epm::RTPField &rRTPField, SmartTruncation pTrunc)
{
   epm::CartesianExpressions   helper(pTrunc);

   int rows = rRTPField.nPh();
   int cols;

   epm::Matrix xComp;
   epm::Matrix yComp;
   epm::Matrix zComp;

   for(int n = 0; n < rRTPField.nR(); n++)
   {
      cols = rRTPField.nTh(n);
      xComp.resize(rows, cols);
      yComp.resize(rows, cols);
      zComp.resize(rows, cols);

      rRTPField.rR().rShell(n).setConstant(0.0);
      rRTPField.rTheta().rShell(n).setConstant(0.0);
      rRTPField.rPhi().rShell(n).setConstant(0.0);

      xComp = -4.0*helper.y(n, 1);
      yComp = helper.z(n, 1) - 12.0*helper.z(n, 2);
      zComp = helper.y(n, 1);

      rRTPField.rR().rShell(n) = helper.xyzToR(xComp, yComp, zComp, n);
      rRTPField.rTheta().rShell(n) = helper.xyzToTh(xComp, yComp, zComp, n);
      rRTPField.rPhi().rShell(n) = helper.xyzToPh(xComp, yComp, zComp, n);
   }
}

/**
 * @brief SpectralSH Toroidal/Poloidal field
 */
void initTorPolSSH(TorPolFieldType &rToPoField)
{
   int nN = rToPoField.trunc()->sim()->rad()->nN();
   int nM;
   int nL = rToPoField.trunc()->local()->spec()->nL();

   for(int l=0; l < nL; ++l)
   {
      nM = rToPoField.trunc()->local()->spec()->nM(l);
      rToPoField.rTor().rLShell(l).setZero(nN, nM);
      rToPoField.rPol().rLShell(l).setZero(nN, nM);
   }

   epm::ArrayI ls = rToPoField.trunc()->local()->spec()->lArray();
   for(int l=rToPoField.tor().minL(); l < nL; ++l)
   {
      rToPoField.rTor().rLShell(l).setConstant(epm::EPMComplex(ls(l),-ls(l)));
      rToPoField.rPol().rLShell(l).setConstant(epm::EPMComplex(ls(l),-ls(l)));

      if(rToPoField.trunc()->local()->spec()->mArray(l)(0) == 0)
      {
         for(int n =0; n< nN; ++n)
         {
            rToPoField.rTor().rLShell(l)(n,0).imag() = 0.0;
            rToPoField.rPol().rLShell(l)(n,0).imag() = 0.0;
         }
      }
   }
}

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
 * @brief Scalar Test starting from RTP
 */
int runRTPScalarTest(SmartTruncation  pTrunc, SSHTransformType &sshTrans)
{
   // Increment number of performed tests
   PERFORMED_TESTS++;

   if(pTrunc->para().id() == 0)
   {
      std::cout << "Scalar precision test starting from RTP" << std::endl;
   }

   int failed = 0;

   // Initialise Error tool
   epm::ErrorComputer  error;

   // Create Scalar test objects
   epm::RTPScalar  rtpScalar(pTrunc);
   epm::RTPScalar  rtpScalar2(pTrunc);
   SpectralScalarType  sshScalar(pTrunc, true);

   // Init RTP scalar values
   initRTPScalar(rtpScalar, pTrunc);

   // Compute transform loop
   for(int i = 0; i < TRANSFORM_STEPS; ++i)
   {
      combineSpectralTransforms(i, sshTrans);
      sshTrans.transformRTP2Spec(sshScalar, rtpScalar);
   }

   for(int i = 0; i < TRANSFORM_STEPS; ++i)
   {
      combineRTPTransforms(i, sshTrans);
      sshTrans.transformSpec2RTP(rtpScalar2, sshScalar);
   }

   for(int n=0; n< pTrunc->local()->rtp()->nR(); ++n)
   {
      failed += error.checkPrecision(rtpScalar.shell(n), rtpScalar2.shell(n));
   }

   epm::Array errMax(2);
   errMax(0) = error.max();
   errMax(1) = error.relMax();

   // Gather max error over all CPUs
   #ifdef EPMDYNAMO_MPI
      // For MPI case the number of CPU is set according to how it's run
      MPI_Allreduce(MPI_IN_PLACE, errMax.data(), 2, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
   #endif //EPMDYNAMO_MPI

   if(pTrunc->para().id() == 0)
   {
      std::cout << "\t Maximum Forward + Backward transforms error: " << errMax(0) << std::endl;
      std::cout << "\t Maximum Forward + Backward transforms relative error: " << errMax(1) << std::endl;
   }

   // Gather total failed comparison from MPI run
   #ifdef EPMDYNAMO_MPI
      // For MPI case the number of CPU is set according to how it's run
      MPI_Allreduce(MPI_IN_PLACE, &failed, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
   #endif //EPMDYNAMO_MPI

   if(pTrunc->para().id() == 0)
   {
      if(failed != 0)
      {
         std::cout << "\t (" << failed << " failed)" << std::endl << std::endl;
      } else
      {
         std::cout << "\t (success)" << std::endl << std::endl;
      }
   }

   return std::min(failed, 1);
}

/**
 * @brief Scalar test starting from SpectralSH
 */
int runSSHScalarTest(SmartTruncation  pTrunc, SSHTransformType &sshTrans)
{
   // Increment number of performed tests
   PERFORMED_TESTS++;

   if(pTrunc->para().id() == 0)
   {
      std::cout << "Scalar precision test starting from SSH" << std::endl;
   }

   int failed = 0;

   // Initialise Error tool
   epm::ErrorComputer  error;

   // Create Scalar test objects
   epm::RTPScalar  rtpScalar(pTrunc);
   SpectralScalarType sshScalar(pTrunc, true);
   SpectralScalarType sshScalar2(pTrunc, true);

   // Init FSHScalar
   initSpectralSHScalar(sshScalar);

   // Compute transform loop
   // Compute transform loop
   for(int i = 0; i < TRANSFORM_STEPS; ++i)
   {
      combineRTPTransforms(i, sshTrans);
      sshTrans.transformSpec2RTP(rtpScalar, sshScalar);
   }

   // Compute transform loop
   for(int i = 0; i < TRANSFORM_STEPS; ++i)
   {
      combineSpectralTransforms(i, sshTrans);
      sshTrans.transformRTP2Spec(sshScalar2, rtpScalar);
   }

   for(int l=0; l < pTrunc->local()->spec()->nL(); ++l)
   {
      failed += error.checkPrecision(sshScalar.lshell(l), sshScalar2.lshell(l));
   }

   epm::Array errMax(2);
   errMax(0) = error.max();
   errMax(1) = error.relMax();

   // Gather max error over all CPUs
   #ifdef EPMDYNAMO_MPI
      // For MPI case the number of CPU is set according to how it's run
      MPI_Allreduce(MPI_IN_PLACE, errMax.data(), 2, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
   #endif //EPMDYNAMO_MPI

   if(pTrunc->para().id() == 0)
   {
      std::cout << "\t Maximum Forward + Backward transforms error: " << errMax(0) << std::endl;
      std::cout << "\t Maximum Forward + Backward transforms relative error: " << errMax(1) << std::endl;
   }

   // Gather total failed comparison from MPI run
   #ifdef EPMDYNAMO_MPI
      // For MPI case the number of CPU is set according to how it's run
      MPI_Allreduce(MPI_IN_PLACE, &failed, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
   #endif //EPMDYNAMO_MPI

   if(pTrunc->para().id() == 0)
   {
      if(failed != 0)
      {
         std::cout << "\t (" << failed << " failed)" << std::endl << std::endl;
      } else
      {
         std::cout << "\t (success)" << std::endl << std::endl;
      }
   }

   return std::min(failed, 1);
}

/**
 * @brief Toroidal/Poloidal test from RTP
 */
int runRTPDiv0Test(SmartTruncation  pTrunc, SSHTransformType &sshTrans)
{
   // Increment number of performed tests
   PERFORMED_TESTS++;

   if(pTrunc->para().id() == 0)
   {
      std::cout << "Toroidal/Poloidal precision test from RTP" << std::endl;
   }

   int failed = 0;

   // Initialise Error tool
   epm::ErrorComputer  error;

   // Create RTP field test objects
   epm::RTPField  rtpField(pTrunc);
   epm::RTPField  rtpField2(pTrunc);

   // Create Toroidal/Poloidal test objects
   TorPolFieldType  topoField(pTrunc);

   //Init RTP field values
   initDiv0RTPField(rtpField, pTrunc);

   // Compute transform loop
   for(int i = 0; i < TRANSFORM_STEPS; ++i)
   {
      combineSpectralTransforms(i, sshTrans);
      sshTrans.transformRTP2TorPol(topoField, rtpField);
   }

   // Compute transform loop
   for(int i = 0; i < TRANSFORM_STEPS; ++i)
   {
      combineRTPTransforms(i, sshTrans);
      sshTrans.transformTorPol2RTP(rtpField2, topoField);
   }

   for(int n=0; n< pTrunc->local()->rtp()->nR(); ++n)
   {
      failed += error.checkPrecision(rtpField.r().shell(n), rtpField2.r().shell(n));
   }

   epm::Array errMax(2);
   errMax(0) = error.max();
   errMax(1) = error.relMax();

   // Gather max error over all CPUs
   #ifdef EPMDYNAMO_MPI
      // For MPI case the number of CPU is set according to how it's run
      MPI_Allreduce(MPI_IN_PLACE, errMax.data(), 2, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
   #endif //EPMDYNAMO_MPI

   if(pTrunc->para().id() == 0)
   {
      std::cout << "\t Maximum Forward + Backward transforms error (R Component) :" << errMax(0) << std::endl;
      std::cout << "\t Maximum Forward + Backward transforms relative error (R Component) :" << errMax(1) << std::endl;
   }
   error.resetErrors();

   for(int n=0; n< pTrunc->local()->rtp()->nR(); ++n)
   {
      failed += error.checkPrecision(rtpField.theta().shell(n), rtpField2.theta().shell(n));
   }

   errMax(0) = error.max();
   errMax(1) = error.relMax();

   // Gather max error over all CPUs
   #ifdef EPMDYNAMO_MPI
      // For MPI case the number of CPU is set according to how it's run
      MPI_Allreduce(MPI_IN_PLACE, errMax.data(), 2, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
   #endif //EPMDYNAMO_MPI

   if(pTrunc->para().id() == 0)
   {
      std::cout << "\t Maximum Forward + Backward transforms error (Theta Component): " << errMax(0) << std::endl;
      std::cout << "\t Maximum Forward + Backward transforms relative error (Theta Component): " << errMax(1) << std::endl;
   }
   error.resetErrors();

   for(int n=0; n< pTrunc->local()->rtp()->nR(); ++n)
   {
      failed += error.checkPrecision(rtpField.phi().shell(n), rtpField2.phi().shell(n));
   }

   errMax(0) = error.max();
   errMax(1) = error.relMax();

   // Gather max error over all CPUs
   #ifdef EPMDYNAMO_MPI
      // For MPI case the number of CPU is set according to how it's run
      MPI_Allreduce(MPI_IN_PLACE, errMax.data(), 2, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
   #endif //EPMDYNAMO_MPI

   if(pTrunc->para().id() == 0)
   {
      std::cout << "\t Maximum Forward + Backward transforms error (Phi Component): " << errMax(0) << std::endl;
      std::cout << "\t Maximum Forward + Backward transforms relative error (Phi Component): " << errMax(1) << std::endl;
   }
   error.resetErrors();

   // Gather total failed comparison from MPI run
   #ifdef EPMDYNAMO_MPI
      // For MPI case the number of CPU is set according to how it's run
      MPI_Allreduce(MPI_IN_PLACE, &failed, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
   #endif //EPMDYNAMO_MPI

   if(pTrunc->para().id() == 0)
   {
      if(failed != 0)
      {
         std::cout << "\t (" << failed << " failed)" << std::endl << std::endl;
      } else
      {
         std::cout << "\t (success)" << std::endl << std::endl;
      }
   }

   return std::min(failed, 1);
}

/**
 * @brief Toroidal/Poloidal test from SpectralSH
 */
int runTorPolTest(SmartTruncation  pTrunc, SSHTransformType &sshTrans)
{
   // Increment number of performed tests
   PERFORMED_TESTS++;

   if(pTrunc->para().id() == 0)
   {
      std::cout << "Toroidal/Poloidal precision test from SpectralSH" << std::endl;
   }

   int failed = 0;

   // Initialise Error tool
   epm::ErrorComputer  error;

   // Create RTP field test objects
   epm::RTPField  rtpField(pTrunc);

   // Create Toroidal/Poloidal test objects
   TorPolFieldType  topoField(pTrunc);
   TorPolFieldType  topoField2(pTrunc);

   //Init SSH field values
   initTorPolSSH(topoField);

   // Compute transform loop
   for(int i = 0; i < TRANSFORM_STEPS; ++i)
   {
      combineRTPTransforms(i, sshTrans);
      sshTrans.transformTorPol2RTP(rtpField, topoField);
   }

   // Compute transform loop
   for(int i = 0; i < TRANSFORM_STEPS; ++i)
   {
      combineSpectralTransforms(i, sshTrans);
      sshTrans.transformRTP2TorPol(topoField2, rtpField);
   }

   for(int l=0; l < pTrunc->local()->spec()->nL(); ++l)
   {
      failed += error.checkPrecision(topoField.tor().lshell(l), topoField2.tor().lshell(l));
   }

   epm::Array errMax(2);
   errMax(0) = error.max();
   errMax(1) = error.relMax();

   // Gather max error over all CPUs
   #ifdef EPMDYNAMO_MPI
      // For MPI case the number of CPU is set according to how it's run
      MPI_Allreduce(MPI_IN_PLACE, errMax.data(), 2, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
   #endif //EPMDYNAMO_MPI

   if(pTrunc->para().id() == 0)
   {
      std::cout << "\t Maximum Backward + Forward transforms error (Toroidal Component): " << errMax(0) << std::endl;
      std::cout << "\t Maximum Backward + Forward transforms relative error (Toroidal Component): " << errMax(1) << std::endl;
   }
   error.resetErrors();

   for(int l=0; l < pTrunc->local()->spec()->nL(); ++l)
   {
      failed += error.checkPrecision(topoField.pol().lshell(l), topoField2.pol().lshell(l));
   }

   errMax(0) = error.max();
   errMax(1) = error.relMax();

   // Gather max error over all CPUs
   #ifdef EPMDYNAMO_MPI
      // For MPI case the number of CPU is set according to how it's run
      MPI_Allreduce(MPI_IN_PLACE, errMax.data(), 2, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
   #endif //EPMDYNAMO_MPI

   if(pTrunc->para().id() == 0)
   {
      std::cout << "\t Maximum Backward + Forward transforms error (Poloidal Component): " << errMax(0) << std::endl;
      std::cout << "\t Maximum Backward + Forward transforms relative error (Poloidal Component): " << errMax(1) << std::endl;
   }
   error.resetErrors();

   // Gather total failed comparison from MPI run
   #ifdef EPMDYNAMO_MPI
      // For MPI case the number of CPU is set according to how it's run
      MPI_Allreduce(MPI_IN_PLACE, &failed, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
   #endif //EPMDYNAMO_MPI

   if(pTrunc->para().id() == 0)
   {
      if(failed != 0)
      {
         std::cout << "\t (" << failed << " failed)" << std::endl << std::endl;
      } else
      {
         std::cout << "\t (success)" << std::endl << std::endl;
      }
   }

   return std::min(failed, 1);
}

/**
 * @brief Gradient test
 */
int runGradientTest(SmartTruncation  pTrunc, SSHTransformType &sshTrans)
{
   // Increment number of performed tests
   PERFORMED_TESTS++;

   if(pTrunc->para().id() == 0)
   {
      std::cout << "Gradient computation check" << std::endl;
   }

   int failed = 0;

   // Initialise Error tool
   epm::ErrorComputer  error;

   // Create Scalar test objects
   epm::RTPScalar  rtpScalar(pTrunc);
   SpectralScalarType sshScalar(pTrunc, true);

   // Create RTP field test objects
   epm::RTPField  rtpField(pTrunc);
   epm::RTPField  rtpField2(pTrunc);

   // Init RTP scalar values
   initRTPGradientScalar(rtpScalar, pTrunc);

   // Compute transform loop
   for(int i = 0; i < TRANSFORM_STEPS; ++i)
   {
      combineSpectralTransforms(i, sshTrans);
      sshTrans.transformRTP2Spec(sshScalar, rtpScalar);
   }

   // Compute transform loop
   for(int i = 0; i < TRANSFORM_STEPS; ++i)
   {
      combineRTPTransforms(i, sshTrans);
      sshTrans.transformSpec2Grad(rtpField2, sshScalar);
   }

   initRTPGradientField(rtpField, pTrunc);

   for(int n=0; n< pTrunc->local()->rtp()->nR(); ++n)
   {
      failed += error.checkPrecision(rtpField.r().shell(n), rtpField2.r().shell(n));
   }

   epm::Array errMax(2);
   errMax(0) = error.max();
   errMax(1) = error.relMax();

   // Gather max error over all CPUs
   #ifdef EPMDYNAMO_MPI
      // For MPI case the number of CPU is set according to how it's run
      MPI_Allreduce(MPI_IN_PLACE, errMax.data(), 2, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
   #endif //EPMDYNAMO_MPI

   if(pTrunc->para().id() == 0)
   {
      std::cout << "\t Maximum Forward + Backward transforms error (R Component) :" << errMax(0) << std::endl;
      std::cout << "\t Maximum Forward + Backward transforms relative error (R Component) :" << errMax(1) << std::endl;
   }
   error.resetErrors();

   for(int n=0; n< pTrunc->local()->rtp()->nR(); ++n)
   {
      failed += error.checkPrecision(rtpField.theta().shell(n), rtpField2.theta().shell(n));
   }

   errMax(0) = error.max();
   errMax(1) = error.relMax();

   // Gather max error over all CPUs
   #ifdef EPMDYNAMO_MPI
      // For MPI case the number of CPU is set according to how it's run
      MPI_Allreduce(MPI_IN_PLACE, errMax.data(), 2, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
   #endif //EPMDYNAMO_MPI

   if(pTrunc->para().id() == 0)
   {
      std::cout << "\t Maximum Forward + Backward transforms error (Theta Component): " << errMax(0) << std::endl;
      std::cout << "\t Maximum Forward + Backward transforms relative error (Theta Component): " << errMax(1) << std::endl;
   }
   error.resetErrors();

   for(int n=0; n< pTrunc->local()->rtp()->nR(); ++n)
   {
      failed += error.checkPrecision(rtpField.phi().shell(n), rtpField2.phi().shell(n));
   }

   errMax(0) = error.max();
   errMax(1) = error.relMax();

   // Gather max error over all CPUs
   #ifdef EPMDYNAMO_MPI
      // For MPI case the number of CPU is set according to how it's run
      MPI_Allreduce(MPI_IN_PLACE, errMax.data(), 2, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
   #endif //EPMDYNAMO_MPI

   if(pTrunc->para().id() == 0)
   {
      std::cout << "\t Maximum Forward + Backward transforms error (Phi Component): " << errMax(0) << std::endl;
      std::cout << "\t Maximum Forward + Backward transforms relative error (Phi Component): " << errMax(1) << std::endl;
   }
   error.resetErrors();

   // Gather total failed comparison from MPI run
   #ifdef EPMDYNAMO_MPI
      // For MPI case the number of CPU is set according to how it's run
      MPI_Allreduce(MPI_IN_PLACE, &failed, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
   #endif //EPMDYNAMO_MPI

   if(pTrunc->para().id() == 0)
   {
      if(failed != 0)
      {
         std::cout << "\t (" << failed << " failed)" << std::endl << std::endl;
      } else
      {
         std::cout << "\t (success)" << std::endl << std::endl;
      }
   }

   return std::min(failed, 1);
}

/**
 * @brief Curl test
 */
int runCurlTest(SmartTruncation  pTrunc, SSHTransformType &sshTrans)
{
   // Increment number of performed tests
   PERFORMED_TESTS++;

   if(pTrunc->para().id() == 0)
   {
      std::cout << "Toroidal/Poloidal curl precision test" << std::endl;
   }

   int failed = 0;

   // Initialise Error tool
   epm::ErrorComputer  error;

   // Create RTP field test objects
   epm::RTPField  rtpField(pTrunc);
   epm::RTPField  rtpField2(pTrunc);

   // Create Toroidal/Poloidal test objects
   TorPolFieldType  topoField(pTrunc);

   //Init RTP field values
   initRTPCurlTorPol(rtpField, pTrunc);

   // Compute transform loop
   for(int i = 0; i < TRANSFORM_STEPS; ++i)
   {
      combineSpectralTransforms(i, sshTrans);
      sshTrans.transformRTP2TorPol(topoField, rtpField);
   }

   // Compute transform loop
   for(int i = 0; i < TRANSFORM_STEPS; ++i)
   {
      combineRTPTransforms(i, sshTrans);
      sshTrans.transformTorPol2Curl(rtpField, topoField);
   }

   initRTPCurlResult(rtpField2, pTrunc);

   for(int n=0; n< pTrunc->local()->rtp()->nR(); ++n)
   {
      failed += error.checkPrecision(rtpField.r().shell(n), rtpField2.r().shell(n));
   }

   epm::Array errMax(2);
   errMax(0) = error.max();
   errMax(1) = error.relMax();

   // Gather max error over all CPUs
   #ifdef EPMDYNAMO_MPI
      // For MPI case the number of CPU is set according to how it's run
      MPI_Allreduce(MPI_IN_PLACE, errMax.data(), 2, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
   #endif //EPMDYNAMO_MPI

   if(pTrunc->para().id() == 0)
   {
      std::cout << "\t Maximum Forward + Backward transforms error (R Component) :" << errMax(0) << std::endl;
      std::cout << "\t Maximum Forward + Backward transforms relative error (R Component) :" << errMax(1) << std::endl;
   }
   error.resetErrors();

   for(int n=0; n< pTrunc->local()->rtp()->nR(); ++n)
   {
      failed += error.checkPrecision(rtpField.theta().shell(n), rtpField2.theta().shell(n));
   }

   errMax(0) = error.max();
   errMax(1) = error.relMax();

   // Gather max error over all CPUs
   #ifdef EPMDYNAMO_MPI
      // For MPI case the number of CPU is set according to how it's run
      MPI_Allreduce(MPI_IN_PLACE, errMax.data(), 2, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
   #endif //EPMDYNAMO_MPI

   if(pTrunc->para().id() == 0)
   {
      std::cout << "\t Maximum Forward + Backward transforms error (Theta Component): " << errMax(0) << std::endl;
      std::cout << "\t Maximum Forward + Backward transforms relative error (Theta Component): " << errMax(1) << std::endl;
   }
   error.resetErrors();

   for(int n=0; n< pTrunc->local()->rtp()->nR(); ++n)
   {
      failed += error.checkPrecision(rtpField.phi().shell(n), rtpField2.phi().shell(n));
   }

   errMax(0) = error.max();
   errMax(1) = error.relMax();

   // Gather max error over all CPUs
   #ifdef EPMDYNAMO_MPI
      // For MPI case the number of CPU is set according to how it's run
      MPI_Allreduce(MPI_IN_PLACE, errMax.data(), 2, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
   #endif //EPMDYNAMO_MPI

   if(pTrunc->para().id() == 0)
   {
      std::cout << "\t Maximum Forward + Backward transforms error (Phi Component): " << errMax(0) << std::endl;
      std::cout << "\t Maximum Forward + Backward transforms relative error (Phi Component): " << errMax(1) << std::endl;
   }
   error.resetErrors();

   // Gather total failed comparison from MPI run
   #ifdef EPMDYNAMO_MPI
      // For MPI case the number of CPU is set according to how it's run
      MPI_Allreduce(MPI_IN_PLACE, &failed, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
   #endif //EPMDYNAMO_MPI

   if(pTrunc->para().id() == 0)
   {
      if(failed != 0)
      {
         std::cout << "\t (" << failed << " failed)" << std::endl << std::endl;
      } else
      {
         std::cout << "\t (success)" << std::endl << std::endl;
      }
   }

   return std::min(failed, 1);
}

/**
 * @brief Nonlinear test
 */
int runNTermsTest(SmartTruncation  pTrunc, SSHTransformType &sshTrans)
{
   // Increment number of performed tests
   PERFORMED_TESTS++;

   if(pTrunc->para().id() == 0)
   {
      std::cout << "Non linear terms projections computation check" << std::endl;
   }

   int failed = 0;

   // Initialise Error tool
   epm::ErrorComputer  error;

   // Create Scalar test objects
   epm::RTPScalar  rtpScalar(pTrunc);
   epm::RTPScalar  rtpScalar2(pTrunc);
   SpectralScalarType sshScalar(pTrunc, true);
   SpectralScalarType sshScalar2(pTrunc, true);

   // Create RTP field test objects
   epm::RTPField  rtpField(pTrunc);
   epm::RTPField  rtpField2(pTrunc);

   // Create Toroidal/Poloidal test objects
   TorPolFieldType  topoField(pTrunc);

   // Init RTP scalar values
   initRTPNTerms(rtpField, pTrunc);

   // Compute transform loop
   for(int i = 0; i < TRANSFORM_STEPS; ++i)
   {
      combineSpectralTransforms(i, sshTrans);
      sshTrans.transformRTP2TorPolNTerms(topoField.rTor(), topoField.rPol(), rtpField);
   }

   initRTPNTermsCurlResult(rtpField2, pTrunc);

   for(int n=0; n< pTrunc->local()->rtp()->nR(); ++n)
   {
      rtpScalar.rShell(n) = rtpField2.r().shell(n) * pTrunc->local()->rtp()->radGrid()(n);
   }

   // Compute transform loop
   for(int i = 0; i < TRANSFORM_STEPS; ++i)
   {
      combineSpectralTransforms(i, sshTrans);
      sshTrans.transformRTP2Spec(sshScalar, rtpScalar);
   }

   initRTPNTermsCurlCurlResult(rtpField2, pTrunc);
   for(int n=0; n< pTrunc->local()->rtp()->nR(); ++n)
   {
      rtpScalar2.rShell(n) = rtpField2.r().shell(n) * pTrunc->local()->rtp()->radGrid()(n);
   }

   // Compute transform loop
   for(int i = 0; i < TRANSFORM_STEPS; ++i)
   {
      combineSpectralTransforms(i, sshTrans);
      sshTrans.transformRTP2Spec(sshScalar2, rtpScalar2);
   }

   epm::ArrayI ls = topoField.trunc()->local()->spec()->lArray();
   for(int l = 0; l < topoField.trunc()->local()->spec()->nL(); ++l)
   {
      if(ls(l) == 0)
      {
         sshScalar.rLShell(l) *= 0.0;
         sshScalar2.rLShell(l) *= 0.0;
         topoField.rTor().rLShell(0).setConstant(0.0);
         topoField.rPol().rLShell(0).setConstant(0.0);
      } else
      {
         sshScalar.rLShell(l) *= 1.0/static_cast<epm::EPMFloat>(ls(l)*(ls(l)+1));
         sshScalar2.rLShell(l) *= 1.0/static_cast<epm::EPMFloat>(ls(l)*(ls(l)+1));
      }
   }

   for(int l=0; l< pTrunc->local()->spec()->nL(); ++l)
   {
      failed += error.checkPrecision(topoField.tor().lshell(l), sshScalar.lshell(l));
   }

   epm::Array errMax(2);
   errMax(0) = error.max();
   errMax(1) = error.relMax();

   // Gather max error over all CPUs
   #ifdef EPMDYNAMO_MPI
      // For MPI case the number of CPU is set according to how it's run
      MPI_Allreduce(MPI_IN_PLACE, errMax.data(), 2, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
   #endif //EPMDYNAMO_MPI

   if(pTrunc->para().id() == 0)
   {
      std::cout << "\t Maximum curl error (Tor Component): " << errMax(0) << std::endl;
      std::cout << "\t Maximum curl relative error (Tor Component): " << errMax(1) << std::endl;
   }
   error.resetErrors();

   for(int l=0; l< pTrunc->local()->spec()->nL(); ++l)
   {
      failed += error.checkPrecision(topoField.pol().lshell(l), sshScalar2.lshell(l));
   }

   errMax(0) = error.max();
   errMax(1) = error.relMax();

   // Gather max error over all CPUs
   #ifdef EPMDYNAMO_MPI
      // For MPI case the number of CPU is set according to how it's run
      MPI_Allreduce(MPI_IN_PLACE, errMax.data(), 2, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
   #endif //EPMDYNAMO_MPI

   if(pTrunc->para().id() == 0)
   {
      std::cout << "\t Maximum curl curl error (Pol Component): " << errMax(0) << std::endl;
      std::cout << "\t Maximum curl curl relative error (Pol Component): " << errMax(1) << std::endl;
   }
   error.resetErrors();

   // Gather total failed comparison from MPI run
   #ifdef EPMDYNAMO_MPI
      // For MPI case the number of CPU is set according to how it's run
      MPI_Allreduce(MPI_IN_PLACE, &failed, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
   #endif //EPMDYNAMO_MPI

   if(pTrunc->para().id() == 0)
   {
      if(failed != 0)
      {
         std::cout << "\t (" << failed << " failed)" << std::endl << std::endl;
      } else
      {
         std::cout << "\t (success)" << std::endl << std::endl;
      }
   }

   return std::min(failed, 1);
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
int runPrecTest()
{
   // Set some truncation values
   int maxN = 20;
   int maxL = 32;
   int maxM = 32;
   int Mp = 1;
   int nCore = 1;

   #ifdef EPMDYNAMO_MPI
      // For MPI case the number of CPU is set according to how it's run
      MPI_Comm_size(MPI_COMM_WORLD, &nCore);
   #endif //EPMDYNAMO_MPI

   // Initialise Spectral truncation information
   SmartTruncation   pTrunc = SimulationType::createTrunc(maxN, maxL, maxM, Mp, nCore);

   if(pTrunc->para().id() == 0)
   {
      std::cout << "Truncation information:" << std::endl;
      std::cout << "\t Radial N: " << pTrunc->sim()->rad()->nR() << std::endl;
      std::cout << "\t Theta N: " << pTrunc->sim()->hoz()->nTh() << std::endl;
      std::cout << "\t Phi N: " << pTrunc->sim()->hoz()->nPh() << std::endl;
      std::cout << "\t Max N: " << pTrunc->sim()->rad()->maxN() << std::endl;
      std::cout << "\t Max L: " << pTrunc->sim()->hoz()->maxL() << std::endl;
      std::cout << "\t Max M: " << pTrunc->sim()->hoz()->maxM() << std::endl;
      std::cout << "\t Mp : " << pTrunc->sim()->hoz()->mp() << std::endl;
      std::cout << std::endl << std::endl;
   }

   // Initialise additional polynomial Basis for initialisation
   epm::AssocLegendreBasis   poly(pTrunc);
   RadialBasisType   radPoly(pTrunc);

   // Initialise the grids of the truncation
   pTrunc->setRadialGrid(radPoly.at(0).pGrid());
   pTrunc->setThetaGrid(poly.at(0).pTheta(), poly.at(0).pCosTheta(), poly.at(0).pSinTheta());

   // Create Transform object
   SSHTransformType    sshTrans(pTrunc);

   // Setup the transform object for parallel versions
   setupTransform(sshTrans);

   int failed = 0;

   // Make sure CPUs are synchronized before start of test
   epm::EPMDYNAMO_SYNCHRONIZE;

   // Run RTP scalar test
   failed += runRTPScalarTest(pTrunc, sshTrans);

   // Make sure CPUs are synchronized before start of test
   epm::EPMDYNAMO_SYNCHRONIZE;

   // Run SSH scalar test
   failed += runSSHScalarTest(pTrunc, sshTrans);

   // Make sure CPUs are synchronized before start of test
   epm::EPMDYNAMO_SYNCHRONIZE;

   // Run RTP div0 test
   failed += runRTPDiv0Test(pTrunc, sshTrans);

   // Make sure CPUs are synchronized before start of test
   epm::EPMDYNAMO_SYNCHRONIZE;

   // Run Toroidal/Poloida test
   failed += runTorPolTest(pTrunc, sshTrans);

   // Make sure CPUs are synchronized before start of test
   epm::EPMDYNAMO_SYNCHRONIZE;

   // Run Gradient test
   failed += runGradientTest(pTrunc, sshTrans);

   // Make sure CPUs are synchronized before start of test
   epm::EPMDYNAMO_SYNCHRONIZE;

   // Run Curl test
   failed += runCurlTest(pTrunc, sshTrans);

   // Make sure CPUs are synchronized before start of test
   epm::EPMDYNAMO_SYNCHRONIZE;

   // Run Non linear terms test
   failed += runNTermsTest(pTrunc, sshTrans);

   // Make sure CPUs are synchronized before gathering data
   epm::EPMDYNAMO_SYNCHRONIZE;

   return failed;
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
      code = runPrecTest();
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
         std::cout << "Failed! (" << code << " out of " << PERFORMED_TESTS << ")" << std::endl;
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
