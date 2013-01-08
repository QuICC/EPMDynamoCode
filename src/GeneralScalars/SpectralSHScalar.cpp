/** \file SpectralSHScalar.cpp
 *  \brief Source of the implementation of spectral space data type
 */

// System includes
//
#include "Config/Parallelisation.h"

// External includes
//

// Class include
//
#include "GeneralScalars/SpectralSHScalar.hpp"

// Project includes
//

namespace EPMDynamo {

   SpectralSHScalar::SpectralSHScalar(const SmartTruncation pTrunc, const bool needsL0)
      : SpectralSHScalarBase(pTrunc, needsL0), mEEven(0.0), mEOdd(0.0), mSpectrumM(pTrunc->sim()->hoz()->nM()), mSpectrumL(pTrunc->sim()->hoz()->nL()), mSpectrumN(pTrunc->sim()->rad()->nN(), pTrunc->sim()->hoz()->nL()), mEFactor(1.0) 
   {
   }

   EPMFloat SpectralSHScalar::l2Norm() const
   {
      // Create temporary storage
      EPMFloat norm = 0.0;
      EPMFloat tmp = 0.0;

      // Get truncation information
      const int l0 = this->minL();
      int nL = this->nL();
      int nM;

      // Loop over harmonic modes
      for(int l = l0; l < nL; ++l)
      {
         nM = this->nM(l);
         for(int m=0; m < nM; ++m)
         {
            // Compute L2 norm            
            tmp = this->lshell(l).col(m).real().dot(this->lshell(l).col(m).real()) + this->lshell(l).col(m).imag().dot(this->lshell(l).col(m).imag());
            norm += tmp;
         }
      }

      // Get the "global" norm for MPI code
      #ifdef EPMDYNAMO_MPI
         MPI_Allreduce(MPI_IN_PLACE, &norm, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      #endif // EPMDYNAMO_MPI

      norm = std::sqrt(norm);

      return norm;
   }

   EPMFloat SpectralSHScalar::modalL2Norm() const
   {
      // Create temporary storage
      EPMFloat norm = 0.0;
      EPMFloat tmp = 0.0;

      // Get truncation information
      const int l0 = this->minL();
      int nL = this->nL();
      int nM;

      EPMFloat tmpMax = 0.0;

      // Loop over harmonic modes
      for(int l = l0; l < nL; ++l)
      {
         nM = this->nM(l);
         for(int m=0; m < nM; ++m)
         {
            // Compute L2 norm            
            tmp = this->lshell(l).col(m).real().dot(this->lshell(l).col(m).real()) + this->lshell(l).col(m).imag().dot(this->lshell(l).col(m).imag());
            norm += std::sqrt(tmp);
            tmpMax = std::max(tmpMax, std::sqrt(tmp));
         }
      }

      // Get the "global" norm for MPI code
      #ifdef EPMDYNAMO_MPI
         MPI_Allreduce(MPI_IN_PLACE, &norm, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      #endif // EPMDYNAMO_MPI

      return norm;
   }

   Array SpectralSHScalar::l2NormByL() const
   {
      // Get truncation information
      const int l0 = this->minL();
      int nL = this->trunc()->sim()->hoz()->nL();
      int nM;
      Array norm(nL);
      norm.setConstant(0.0);
      ArrayI ls = this->trunc()->local()->spec()->lArray();

      EPMFloat tmp = 0.0;
      EPMFloat tmpMax = 0.0;

      // Loop over harmonic modes
      for(int i = l0; i < this->nL(); ++i)
      {
         nM = this->nM(i);
         tmpMax = 0.0;
         for(int m=0; m < nM; ++m)
         {
            // Compute L2 norm            
            tmp = this->lshell(i).col(m).real().dot(this->lshell(i).col(m).real()) + this->lshell(i).col(m).imag().dot(this->lshell(i).col(m).imag());
            tmpMax = std::max(tmpMax, std::sqrt(tmp));
         }
         norm(ls(i)) = tmpMax;

      }

      // Get the "global" norm for MPI code
      #ifdef EPMDYNAMO_MPI
         MPI_Allreduce(MPI_IN_PLACE, &norm, norm.size(), MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
      #endif // EPMDYNAMO_MPI

      return norm;
   }

   Array SpectralSHScalar::energy() const
   {
      Array energy(1);

      energy(0) = this->mSpectrumL.sum();

      return energy;
   }

   EPMFloat SpectralSHScalar::evenEnergy() const
   {
      return this->mEEven;
   }

   EPMFloat SpectralSHScalar::oddEnergy() const
   {
      return this->mEOdd;
   }

   void SpectralSHScalar::setEnergyScale(const EPMFloat& eFactor)
   {
      // Se the energy normalisation
      this->mEFactor = eFactor;
   }

   void SpectralSHScalar::scaleSpectra()
   {
      // Check if scaling is required
      if(this->mEFactor != 1.0)
      {
         // Rescale the M spectrum
         this->mSpectrumM *= this->mEFactor;

         // Rescale the L spectrum
         this->mSpectrumL *= this->mEFactor;

         // Rescale the N spectrum
         this->mSpectrumN *= this->mEFactor;
      }
   }

   void SpectralSHScalar::initialiseZeros()
   {
      for(int l=0; l < this->nL(); ++l)
      {
         this->rLShell(l).setConstant(EPMComplex(0.0,0.0));
      }
   }

   void SpectralSHScalar::checkInitialisation()
   {
      int l0 = this->minL();
      int nL = this->nL();

      for(int l = l0; l < nL; ++l)
      {
         if(this->trunc()->local()->spec()->mArray(l)(0) == 0)
         {
            this->rLShell(l).col(0).imag().setConstant(0.0);
         }
      }
   }

}
