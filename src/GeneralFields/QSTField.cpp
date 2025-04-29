/** \file QSTField.cpp
 *  \brief Implementation of the QST field decomposition
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "GeneralFields/QSTField.hpp"

// Project includes
//

namespace EPMDynamo {

   QSTField::QSTField(SmartTruncation pTrunc, const bool needsL0)
      :SpectralFieldBase(pTrunc), mQ(pTrunc, needsL0), mS(pTrunc, needsL0), mT(pTrunc, needsL0)
   {
   }

   void QSTField::rescale(const EPMFloat scale)
   {
      // Rescale Q component
      this->mQ.rescale(scale);

      // Rescale S component
      this->mS.rescale(scale);

      // Rescale T component
      this->mT.rescale(scale);
   }

   void QSTField::computeQSpectra(const QSTField::RadialBasisType &radBasis)
   {
      this->rQ().rSpectrumL().setConstant(0.0);
      this->rQ().rSpectrumM().setConstant(0.0);

      int l0 = this->q().minL();
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
                  tmpEnergy += radBasis.at(l).qsEWeights()(k,n) * (this->q().lshell(l)(n,m).real()*this->q().lshell(l)(k,m).real()+this->q().lshell(l)(n,m).imag()*this->q().lshell(l)(k,m).imag());
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
            this->rQ().rSpectrumM()(ms(m)) += tmpSpectrum(ms(m));
         }

         // Compute energy in L spectrum
         for(int i = 0; i < ms.size(); ++i)
         {
            this->rQ().rSpectrumL()(l_) += tmpSpectrum(ms(i));
         }
      }

      // Get the "global" spectra for MPI code
#ifdef EPMDYNAMO_MPI
      MPI_Allreduce(MPI_IN_PLACE, this->rQ().rSpectrumL().data(), this->q().spectrumL().size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(MPI_IN_PLACE, this->rQ().rSpectrumM().data(), this->q().spectrumM().size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif // EPMDYNAMO_MPI
   }

   void QSTField::computeSSpectra(const QSTField::RadialBasisType &radBasis)
   {
      this->rS().rSpectrumL().setConstant(0.0);
      this->rS().rSpectrumM().setConstant(0.0);

      int l0 = this->s().minL();
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
                  tmpEnergy += radBasis.at(l).qsEWeights()(k,n) * (this->s().lshell(l)(n,m).real()*this->s().lshell(l)(k,m).real()+this->s().lshell(l)(n,m).imag()*this->s().lshell(l)(k,m).imag());
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
            this->rS().rSpectrumM()(ms(m)) += tmpSpectrum(ms(m));
         }

         // Compute energy in L spectrum
         for(int i = 0; i < ms.size(); ++i)
         {
            this->rS().rSpectrumL()(l_) += tmpSpectrum(ms(i));
         }
      }

      // Get the "global" spectra for MPI code
#ifdef EPMDYNAMO_MPI
      MPI_Allreduce(MPI_IN_PLACE, this->rS().rSpectrumL().data(), this->s().spectrumL().size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(MPI_IN_PLACE, this->rS().rSpectrumM().data(), this->s().spectrumM().size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif // EPMDYNAMO_MPI
   }

   void QSTField::computeTSpectra(const QSTField::RadialBasisType &radBasis)
   {
      this->rT().rSpectrumL().setConstant(0.0);
      this->rT().rSpectrumM().setConstant(0.0);

      int l0 = this->t().minL();
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
                  tmpEnergy += radBasis.at(l).eWeights()(k,n) * (this->t().lshell(l)(n,m).real()*this->t().lshell(l)(k,m).real()+this->t().lshell(l)(n,m).imag()*this->t().lshell(l)(k,m).imag());
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
            this->rT().rSpectrumM()(ms(m)) += tmpSpectrum(ms(m));
         }

         // Compute energy in L spectrum
         for(int i = 0; i < ms.size(); ++i)
         {
            this->rT().rSpectrumL()(l_) += tmpSpectrum(ms(i));
         }
      }

      // Get the "global" spectra for MPI code
#ifdef EPMDYNAMO_MPI
      MPI_Allreduce(MPI_IN_PLACE, this->rT().rSpectrumL().data(), this->t().spectrumL().size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(MPI_IN_PLACE, this->rT().rSpectrumM().data(), this->t().spectrumM().size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif // EPMDYNAMO_MPI
   }

   void QSTField::setEnergyScale(const EPMFloat& eFactor)
   {
      // Set normalisation for the Q component
      this->rQ().setEnergyScale(eFactor);

      // Set normalisation for the S component
      this->rS().setEnergyScale(eFactor);

      // Set normalisation for the T component
      this->rT().setEnergyScale(eFactor);
   }

   void QSTField::initialiseZeros()
   {
      // initialise Q component to zeros
      this->rQ().initialiseZeros();

      // initialise S component to zeros
      this->rS().initialiseZeros();

      // initialise T component to zeros
      this->rT().initialiseZeros();
   }
}
