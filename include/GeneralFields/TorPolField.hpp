/** \file TorPolField.hpp
 *  \brief Implementation of the Toroidal/Poloidal vector field decomposition
 */

#ifndef TORPOLFIELD_HPP
#define TORPOLFIELD_HPP

// System includes
//
#include "Config/Parallelisation.h"

// External includes
//

// Project includes
//
#include "General/EPMTypedefs.hpp"
#include "General/MathConstants.hpp"
#include "GeneralFields/SpectralFieldBase.hpp"
#include "Simulations/Traits/SimulationTraits.hpp"
#include "Domain/Truncation.hpp"

namespace EPMDynamo {

   /**
    * \brief Implementation of the Toroidal/Poloidal vector field decomposition
    *
    * \tparam TSimType Type of the simulation
    */
   template <typename TSimType> class TorPolField: public SpectralFieldBase
   {
      /// Typedef from Simulation trait to local truncation type
      typedef typename TSimType::ScalarType    ScalarType;

      /// Typedef from Simulation trait to local truncation type
      typedef typename TSimType::RadialBasisType    RadialBasisType;

      public:
         /**
         * @brief Constructs the field using the given spectral truncation
         *
         * @param pTrunc Smart pointer to truncation information
         * @param needsL0  Needs to consider L=0 mode ?
         */
         TorPolField(SmartTruncation pTrunc, const bool needsL0 = false);

         /**
         * @brief Simple empty destructor
         */
         virtual ~TorPolField() {};

         /**
          * @brief Get the Toroidal component of the field
          */
         const ScalarType&  tor() const;

         /**
          * @brief Set the Toroidal component of the field
          */
         ScalarType&  rTor();

         /**
          * @brief Get the Poloidal component of the field
          */
         const ScalarType&  pol() const;

         /**
          * @brief Set the Poloidal component of the field
          */
         ScalarType&  rPol();

         /**
          * @brief Compute the toroidal energy spectra
          *
          * @param radBasis Radial basis 
          */
         void computeTorSpectra(const RadialBasisType &radBasis);

         /**
          * @brief Compute the poloidal energy spectra
          *
          * @param radBasis Radial basis 
          */
         void computePolSpectra(const RadialBasisType &radBasis);

         /**
          * @brief Get the energy in the toroidal component
          */
         Array  torEnergy() const;

         /**
          * @brief Get the energy in the poloidal component
          */
         Array  polEnergy() const;

         /**
          * @brief Get the energy spectrum in the toroidal component
          */
         const Array&  torSpectrumL() const;

         /**
          * @brief Get the energy spectrum in the poloidal component
          */
         const Array&  polSpectrumL() const;

         /**
          * @brief Get the energy spectrum in the toroidal component
          */
         const Array&  torSpectrumM() const;

         /**
          * @brief Get the energy spectrum in the poloidal component
          */
         const Array&  polSpectrumM() const;
         
      protected:

      private:
         /**
          * @brief Storage for the Toroidal component
          */
         ScalarType   mTor;

         /**
          * @brief Storage for the Poloidal component
          */
         ScalarType   mPol;
   };

   template <typename TSimType> TorPolField<TSimType>::TorPolField(SmartTruncation pTrunc, const bool needsL0)
      : SpectralFieldBase(pTrunc), mTor(pTrunc, needsL0), mPol(pTrunc, needsL0)
   {
   }

   template <typename TSimType> inline const typename TorPolField<TSimType>::ScalarType& TorPolField<TSimType>::tor() const
   {
      return this->mTor;
   }

   template <typename TSimType> inline  typename TorPolField<TSimType>::ScalarType& TorPolField<TSimType>::rTor()
   {
      return this->mTor;
   }

   template <typename TSimType> inline const  typename TorPolField<TSimType>::ScalarType& TorPolField<TSimType>::pol() const
   {
      return this->mPol;
   }

   template <typename TSimType> inline  typename TorPolField<TSimType>::ScalarType& TorPolField<TSimType>::rPol()
   {
      return this->mPol;
   }

   template <typename TSimType> void TorPolField<TSimType>::computeTorSpectra(const typename TorPolField<TSimType>::RadialBasisType &radBasis)
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
         this->rTor().rSpectrumL()(l_) += tmpSpectrum.sum();
      }

      // Get the "global" spectra for MPI code
      #ifdef EPMDYNAMO_MPI
         MPI_Allreduce(MPI_IN_PLACE, this->rTor().rSpectrumL().data(), this->tor().spectrumL().size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
         MPI_Allreduce(MPI_IN_PLACE, this->rTor().rSpectrumM().data(), this->tor().spectrumM().size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      #endif // EPMDYNAMO_MPI
   }

   template <typename TSimType> void TorPolField<TSimType>::computePolSpectra(const typename TorPolField<TSimType>::RadialBasisType &radBasis)
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
         this->rPol().rSpectrumL()(l_) += tmpSpectrum.sum();
      }

      // Get the "global" spectra for MPI code
      #ifdef EPMDYNAMO_MPI
         MPI_Allreduce(MPI_IN_PLACE, this->rPol().rSpectrumL().data(), this->pol().spectrumL().size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
         MPI_Allreduce(MPI_IN_PLACE, this->rPol().rSpectrumM().data(), this->pol().spectrumM().size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      #endif // EPMDYNAMO_MPI
   }
}

#endif // TORPOLFIELD_HPP
