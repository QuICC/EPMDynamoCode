/** \file TorPolField.cpp
 *  \brief Implementation of the Toroidal/Poloidal vector field decomposition
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "GeneralFields/TorPolField.hpp"

// Project includes
//

namespace EPMDynamo {

   TorPolField::TorPolField(SmartTruncation pTrunc, const bool needsL0)
      : SpectralFieldBase(pTrunc), mTor(pTrunc, needsL0), mPol(pTrunc, needsL0)
   {
   }

   void TorPolField::rescale(const EPMFloat scale)
   {
      // Rescale toroidal component
      this->mTor.rescale(scale);

      // Rescale poloidal component
      this->mPol.rescale(scale);
   }

   void TorPolField::computeTorSpectra(const TorPolField::RadialBasisType &radBasis)
   {
      this->rTor().rSpectrumL().setConstant(0.0);
      this->rTor().rSpectrumM().setConstant(0.0);

      int l0 = this->tor().minL();
      int nL = this->nL();
      int nN = this->nN();
      ArrayI   ls = this->trunc()->local()->spec()->lArray();
      ArrayI   ms;
      Array    tmpSpectrum(this->trunc()->sim()->hoz()->nM());
      tmpSpectrum.setConstant(0.0);

      EPMFloat tmpEnergy;
      EPMFloat shWeight;
      EPMFloat shFactor;
      EPMFloat lfactor;
      int l_;
      for(int l = l0; l < nL; ++l)
      {
         l_ = ls(l);
         ms = this->trunc()->local()->spec()->mArray(l);
         shWeight = 4.0 * MathConstants::PI / static_cast<EPMFloat>(2*l_+1);
         lfactor = static_cast<EPMFloat>(l_*(l_+1));

         for(int m = 0; m < this->nM(l); ++m)
         {
            tmpEnergy = 0.0;
            for(int n = 0; n < nN; ++n)
            {
               for(int k = 0; k < nN; ++k)
               {
                  tmpEnergy += radBasis.at(l).eWeights()(k,n) * (this->tor().lshell(l)(n,m).real()*this->tor().lshell(l)(k,m).real()+this->tor().lshell(l)(n,m).imag()*this->tor().lshell(l)(k,m).imag());
               }
            }

            if(ms(m) != 0)
            {
               shFactor = 4.0*shWeight*lfactor;
            } else
            {
               shFactor = shWeight*lfactor;
            }

            tmpSpectrum(ms(m)) = shFactor*tmpEnergy;
            this->rTor().rSpectrumM()(ms(m)) += tmpSpectrum(ms(m));
         }

         // Compute energy in L spectrum
         for(int i = 0; i < ms.size(); ++i)
         {
            this->rTor().rSpectrumL()(l_) += tmpSpectrum(ms(i));
         }
      }

      // Use energy factor to scale energy correctly
      this->rTor().scaleSpectra();

      // Get the "global" spectra for MPI code
      #ifdef EPMDYNAMO_MPI
         MPI_Allreduce(MPI_IN_PLACE, this->rTor().rSpectrumL().data(), this->tor().spectrumL().size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
         MPI_Allreduce(MPI_IN_PLACE, this->rTor().rSpectrumM().data(), this->tor().spectrumM().size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      #endif // EPMDYNAMO_MPI
   }

   void TorPolField::computePolSpectra(const TorPolField::RadialBasisType &radBasis)
   {
      this->rPol().rSpectrumL().setConstant(0.0);
      this->rPol().rSpectrumM().setConstant(0.0);

      int l0 = this->pol().minL();
      int nL = this->nL();
      int nN = this->nN();
      ArrayI   ls = this->trunc()->local()->spec()->lArray();
      ArrayI   ms;
      Array    tmpSpectrum(this->trunc()->sim()->hoz()->nM());
      tmpSpectrum.setConstant(0.0);

      EPMFloat tmpEnergy;
      EPMFloat shWeight;
      EPMFloat shFactor;
      int l_;

      for(int l = l0; l < nL; ++l)
      {
         l_ = ls(l);
         ms = this->trunc()->local()->spec()->mArray(l);
         shWeight = 4.0 * MathConstants::PI / static_cast<EPMFloat>(2*l_+1);

         for(int m = 0; m < this->nM(l); ++m)
         {
            tmpEnergy = 0.0;
            for(int n = 0; n < nN; ++n)
            {
               for(int k = 0; k < nN; ++k)
               {
                  tmpEnergy += radBasis.at(l).polEWeights()(k,n) * (this->pol().lshell(l)(n,m).real()*this->pol().lshell(l)(k,m).real()+this->pol().lshell(l)(n,m).imag()*this->pol().lshell(l)(k,m).imag());
               }
            }

            if(ms(m) != 0)
            {
               shFactor = 4.0*shWeight;
            } else
            {
               shFactor = shWeight;
            }

            tmpSpectrum(ms(m)) = shFactor*tmpEnergy;

            this->rPol().rSpectrumM()(ms(m)) += tmpSpectrum(ms(m));
         }

         // Compute energy in L spectrum
         for(int i = 0; i < ms.size(); ++i)
         {
            this->rPol().rSpectrumL()(l_) += tmpSpectrum(ms(i));
         }
      }

      // Use energy factor to scale energy correctly
      this->rPol().scaleSpectra();

      // Get the "global" spectra for MPI code
      #ifdef EPMDYNAMO_MPI
         MPI_Allreduce(MPI_IN_PLACE, this->rPol().rSpectrumL().data(), this->pol().spectrumL().size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
         MPI_Allreduce(MPI_IN_PLACE, this->rPol().rSpectrumM().data(), this->pol().spectrumM().size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      #endif // EPMDYNAMO_MPI
   }

   void TorPolField::setEnergyScale(const EPMFloat& eFactor)
   {
      // Set normalisation for the toroidal component
      this->rTor().setEnergyScale(eFactor);

      // Set normalisation for the poloidal component
      this->rPol().setEnergyScale(eFactor);
   }

   void TorPolField::initialiseZeros()
   {
      // initialise Toroidal scalar to zero
      this->rTor().initialiseZeros();

      // initialise Poloidal scalar to zero
      this->rPol().initialiseZeros();
   }
}
