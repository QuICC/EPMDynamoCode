/** \file SpectralSHScalar.hpp
 *  \brief Implementation of a spectral space scalar field data type
 */

#ifndef SPECTRALSHSCALAR_HPP
#define SPECTRALSHSCALAR_HPP

// System includes
//
#include <vector>

// External includes
//

// Project includes
//
#include "General/EPMTypedefs.hpp"
#include "General/MathConstants.hpp"
#include "Domain/Truncation.hpp"
#include "GeneralScalars/SpectralSHScalarBase.hpp"
#include "Polynomials/RadialBasis.hpp"

namespace EPMDynamo {

   /**
    * @brief Basis class for representing a scalar field in full spectral space
    *
    * This datatype is used for storing the full spectral expansion of the scalar field.
    * It will also be the building block of the spectral vector fields.
    */
   class SpectralSHScalar: public SpectralSHScalarBase
   {
      public:
         /**
         * @brief Construct the scalar from the given spectral truncation information
         *
         * @param pTrunc    Truncation information
         * @param needsL0   Should computation consider L=0 modes
         */
         SpectralSHScalar(const SmartTruncation pTrunc, const bool needsL0);

         /**
         * @brief Simple empty destructor
         */
         virtual ~SpectralSHScalar() {};

         /**
          * @brief Compute the L2 norm of the spectral values
          */
         EPMFloat l2Norm() const;

         /**
          * @brief Compute an L2 norm of the spectral values for each degree
          *
          * This is an attempt to obtain a norm usable in the "spectral" CFL
          * condition
          */
         Array l2NormByL() const;

         /**
          * @brief Compute the L2 norm for each mode separately
          */
         EPMFloat modalL2Norm() const;

         /**
          * @brief Get the energy of the scalar
          */
         Array energy() const;

         /**
          * @brief Get the energy in the even degrees of the scalar
          */
         EPMFloat evenEnergy() const;

         /**
          * @brief Get the energy in the odd degrees of the scalar
          */
         EPMFloat oddEnergy() const;

         /**
          * @brief Set the (full) energy spectrum on M
          */
         Array& rSpectrumM();

         /**
          * @brief Set the (full) energy spectrum on L
          */
         Array& rSpectrumL();

         /**
          * @brief Set the (full) energy spectrum on N
          */
         Matrix& rSpectrumN();

         /**
          * @brief Get the (full) energy spectrum on M
          */
         const Array& spectrumM() const;

         /**
          * @brief Get the (full) energy spectrum on L
          */
         const Array& spectrumL() const;

         /**
          * @brief Get the (full) energy spectrum on N
          */
         const Matrix& spectrumN() const;

         /**
          * @brief Compute energy spectra
          *
          * @param radBasis Radial basis
          */
         template <typename TPolynomial> void computeSpectra(const RadialBasis<TPolynomial> &radBasis);

         /**
          * @brief Define the non-dimensionlisation factors
          *
          * @param eFactor Energy multiplicative factor
          */
         void setEnergyScale(const EPMFloat&  eFactor);

         /**
          * @brief Scale energy spectra
          */
         void scaleSpectra();

         /**
          * @brief Initialise values to zero
          */
         void initialiseZeros();

         /**
          * @brief Check field set during initialisation
          */
         void checkInitialisation();

      protected:

      private:
         /**
          * @brief Storage for the energy spectrum in M
          */
         Array mSpectrumM;

         /**
          * @brief Storage for the energy spectrum in L
          */
         Array mSpectrumL;

         /**
          * @brief Storage for the energy spectrum in N
          */
         Matrix mSpectrumN;

         /**
          * @brief Energy non-dimensionalisation factor
          */
         EPMFloat mEFactor;
   };

   inline const Array& SpectralSHScalar::spectrumL() const
   {
      return this->mSpectrumL;
   }

   inline const Array& SpectralSHScalar::spectrumM() const
   {
      return this->mSpectrumM;
   }

   inline const Matrix& SpectralSHScalar::spectrumN() const
   {
      return this->mSpectrumN;
   }

   inline Array& SpectralSHScalar::rSpectrumL()
   {
      return this->mSpectrumL;
   }

   inline Array& SpectralSHScalar::rSpectrumM()
   {
      return this->mSpectrumM;
   }

   inline Matrix& SpectralSHScalar::rSpectrumN()
   {
      return this->mSpectrumN;
   }

   template <typename TPolynomial> void SpectralSHScalar::computeSpectra(const RadialBasis<TPolynomial> &radBasis)
   {
      this->mSpectrumN.setConstant(0.0);
      this->mSpectrumL.setConstant(0.0);
      this->mSpectrumM.setConstant(0.0);

      int l0 = this->minL();
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
                  nSpec(n) += 2.0*radBasis.at(l).eWeights()(k,n) * (this->lshell(l)(n,m).real()*this->lshell(l)(k,m).real()+this->lshell(l)(n,m).imag()*this->lshell(l)(k,m).imag());
               }

               nSpec(n) += radBasis.at(l).eWeights()(n,n) * (this->lshell(l)(n,m).real()*this->lshell(l)(n,m).real()+this->lshell(l)(n,m).imag()*this->lshell(l)(n,m).imag());

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
            this->rSpectrumM()(ms(m)) += tmpSpectrum(ms(m));

            this->rSpectrumN().col(l_) += shFactor*nSpec;
         }

         // Compute energy in L spectrum
         for(int i = 0; i < ms.size(); ++i)
         {
            this->rSpectrumL()(l_) += tmpSpectrum(ms(i));
         }
      }

      // Use energy factor to scale energy correctly
      this->scaleSpectra();

      // Get the "global" spectra for MPI code
      #ifdef EPMDYNAMO_MPI
         MPI_Allreduce(MPI_IN_PLACE, this->mSpectrumN.data(), this->mSpectrumN.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
         MPI_Allreduce(MPI_IN_PLACE, this->mSpectrumL.data(), this->mSpectrumL.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
         MPI_Allreduce(MPI_IN_PLACE, this->mSpectrumM.data(), this->mSpectrumM.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      #endif // EPMDYNAMO_MPI
   }

}

#endif // SPECTRALSHSCALAR_HPP
