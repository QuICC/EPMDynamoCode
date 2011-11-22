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
      this->rTor().rSpectrumN().setConstant(0.0);
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
      Array nSpec = Array::Zero(nN);
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
            nSpec.setConstant(0.0);
            for(int n = 0; n < nN; ++n)
            {
               for(int k = 0; k < n; ++k)
               {
                  nSpec(n) += 2.0*radBasis.at(l).eWeights()(k,n) * (this->tor().lshell(l)(n,m).real()*this->tor().lshell(l)(k,m).real()+this->tor().lshell(l)(n,m).imag()*this->tor().lshell(l)(k,m).imag());
               }

               nSpec(n) += radBasis.at(l).eWeights()(n,n) * (this->tor().lshell(l)(n,m).real()*this->tor().lshell(l)(n,m).real()+this->tor().lshell(l)(n,m).imag()*this->tor().lshell(l)(n,m).imag());

               tmpEnergy += nSpec(n);
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

            this->rTor().rSpectrumN().col(l_) += shFactor*nSpec;
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
         MPI_Allreduce(MPI_IN_PLACE, this->rTor().rSpectrumN().data(), this->tor().spectrumN().size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
         MPI_Allreduce(MPI_IN_PLACE, this->rTor().rSpectrumL().data(), this->tor().spectrumL().size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
         MPI_Allreduce(MPI_IN_PLACE, this->rTor().rSpectrumM().data(), this->tor().spectrumM().size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      #endif // EPMDYNAMO_MPI
   }

   void TorPolField::computePolSpectra(const TorPolField::RadialBasisType &radBasis)
   {
      this->rPol().rSpectrumN().setConstant(0.0);
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
      Array nSpec = Array::Zero(nN);
      int l_;

      for(int l = l0; l < nL; ++l)
      {
         l_ = ls(l);
         ms = this->trunc()->local()->spec()->mArray(l);
         shWeight = 4.0 * MathConstants::PI / static_cast<EPMFloat>(2*l_+1);

         for(int m = 0; m < this->nM(l); ++m)
         {
            tmpEnergy = 0.0;
            nSpec.setConstant(0.0);
            for(int n = 0; n < nN; ++n)
            {
               for(int k = 0; k < n; ++k)
               {
                  nSpec(n) += 2.0*radBasis.at(l).polEWeights()(k,n) * (this->pol().lshell(l)(n,m).real()*this->pol().lshell(l)(k,m).real()+this->pol().lshell(l)(n,m).imag()*this->pol().lshell(l)(k,m).imag());
               }

               nSpec(n) += radBasis.at(l).polEWeights()(n,n) * (this->pol().lshell(l)(n,m).real()*this->pol().lshell(l)(n,m).real()+this->pol().lshell(l)(n,m).imag()*this->pol().lshell(l)(n,m).imag());

               tmpEnergy += nSpec(n);
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

            this->rPol().rSpectrumN().col(l_) += shFactor*nSpec;
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
         MPI_Allreduce(MPI_IN_PLACE, this->rPol().rSpectrumN().data(), this->pol().spectrumN().size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
         MPI_Allreduce(MPI_IN_PLACE, this->rPol().rSpectrumL().data(), this->pol().spectrumL().size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
         MPI_Allreduce(MPI_IN_PLACE, this->rPol().rSpectrumM().data(), this->pol().spectrumM().size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      #endif // EPMDYNAMO_MPI
   }

   Array TorPolField::computeXYZSolidProjection(const TorPolField::RadialBasisType &radBasis)
   {
      // Create storage for amplitudes
      Array amplitude = Array::Zero(4);

      // Get the stored ls
      ArrayI ls = this->trunc()->local()->spec()->lArray();
      ArrayI ms;
      int l_ = -1;
      int m_ = -1;

      // Setup the x axis solid body rotation (has only a real component)
      Array xAxis = Array::Zero(this->trunc()->sim()->rad()->nN());
      xAxis(0) = -0.443113;
      EPMFloat xProj = 0.0;

      // Setup the y axis solid body rotation (has only an imaginary component)
      Array yAxis = Array::Zero(this->trunc()->sim()->rad()->nN());
      yAxis(0) = 0.443113;
      EPMFloat yProj = 0.0;

      // Setup the y axis solid body rotation (has only a real component)
      Array zAxis = Array::Zero(this->trunc()->sim()->rad()->nN());
      zAxis(0) = 0.886227;
      EPMFloat zProj = 0.0;

      // Loop over all stored l's
      for(int l = 0; l < ls.size(); l++)
      {
         l_ = ls(l);
         if(l_ == 1)
         {
            ms = this->trunc()->local()->spec()->mArray(l);
            for(int m =0; m < ms.size(); m++)
            {
               m_ = ms(m);
               if(m_ == 1)
               {
                  int nN = this->trunc()->sim()->rad()->nN();

                  for(int n = 0; n < nN; ++n)
                  {
                     for(int k = 0; k < nN; ++k)
                     {
                        xProj += radBasis.at(l).eWeights()(k,n) * (this->tor().lshell(l)(n,m).real()*xAxis(k));
                        yProj += radBasis.at(l).eWeights()(k,n) * ( this->tor().lshell(l)(n,m).imag()*yAxis(k));
                     }
                  }

                  xProj /= xAxis.dot(radBasis.at(l).eWeights().col(0));
                  yProj /= yAxis.dot(radBasis.at(l).eWeights().col(0));

                  amplitude(0) += xProj;
                  amplitude(1) += yProj;
               }

               if(m_ == 0)
               {
                  int nN = this->trunc()->sim()->rad()->nN();

                  for(int n = 0; n < nN; ++n)
                  {
                     for(int k = 0; k < nN; ++k)
                     {
                        zProj += radBasis.at(l).eWeights()(k,n) * (this->tor().lshell(l)(n,m).real()*zAxis(k));
                     }
                  }

                  zProj /= zAxis.dot(radBasis.at(l).eWeights().col(0));

                  amplitude(2) += zProj;
               }
            }
         }
      }

      #ifdef EPMDYNAMO_MPI
         MPI_Allreduce(MPI_IN_PLACE, amplitude.data(), amplitude.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      #endif // EPMDYNAMO_MPI

      amplitude(3) = std::sqrt(amplitude(0)*amplitude(0) + amplitude(1)*amplitude(1)+ amplitude(2)*amplitude(2));

      return amplitude;
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

   void TorPolField::checkInitialisation()
   {
      // check initialisation of Toroidal scalar
      this->rTor().checkInitialisation();

      // check initialisation of Poloidal scalar
      this->rPol().checkInitialisation();

      // Get the lowest harmonic degree
      int l0 = this->trunc()->local()->spec()->lArray()(0);

      // Set l=0 component to zero
      if(l0 == 0)
      {
         // Toroidal component
         this->rTor().rLShell(0).setConstant(EPMComplex(0.0,0.0));

         // Poloidal component
         this->rPol().rLShell(0).setConstant(EPMComplex(0.0,0.0));
      }
   }
}
