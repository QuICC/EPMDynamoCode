/** \file QSTField.hpp
 *  \brief Implementation of the QST field decomposition
 */

#ifndef QSTFIELD_HPP
#define QSTFIELD_HPP

// System includes
//
#include "Config/Parallelisation.h"

// External includes
//

// Project includes
//
#include "General/EPMTypedefs.hpp"
#include "Domain/Truncation.hpp"
#include "GeneralFields/SpectralFieldBase.hpp"
#include "Simulations/Traits/SimulationTraits.hpp"
#include "General/MathConstants.hpp"

namespace EPMDynamo {

   /**
    * \brief Implementation of the QST field decomposition
    *
    * \tparam TSimType Type of the simulation
    */
   template <typename TSimType> class QSTField: public SpectralFieldBase
   {
      /// Typedef from Simulation trait to local truncation type
      typedef typename TSimType::ScalarType    ScalarType;

      /// Typedef from Simulation trait to local truncation type
      typedef typename TSimType::RadialBasisType    RadialBasisType;

      public:
         /**
         * @brief Constructs the field using the given spectral truncation information
         *
         * @param pTrunc Truncation information pointer
         */
         QSTField(SmartTruncation pTrunc);

         /**
         * @brief Simple empty constructor
         */
         virtual ~QSTField() {};

         /**
          * @brief Get Q component of the field
          */
         const ScalarType&  q() const;

         /**
          * @brief Set Q component of the field
          */
         ScalarType&  rQ();

         /**
          * @brief Get S component of the field
          */
         const ScalarType&  s() const;

         /**
          * @brief Set S component of the field
          */
         ScalarType&  rS();

         /**
          * @brief Get T component of the field
          */
         const ScalarType&  t() const;

         /**
          * @brief Set T component of the field
          */
         ScalarType&  rT();

         /**
          * @brief Compute the energy spectra of the Q QST Component
          *
          * \param radBasis Radial basis
          */
         void  computeQSpectra(const RadialBasisType &radBasis);

         /**
          * @brief Compute the energy spectra of the S QST Component
          *
          * \param radBasis Radial basis
          */
         void  computeSSpectra(const RadialBasisType &radBasis);

         /**
          * @brief Compute the energy spectra of the T QST Component
          *
          * \param radBasis Radial basis
          */
         void  computeTSpectra(const RadialBasisType &radBasis);
         
      protected:

      private:
         /**
          * @brief Q component of the QST field
          */
         ScalarType   mQ;

         /**
          * @brief S component of the QST field
          */
         ScalarType   mS;

         /**
          * @brief T component of the QST field
          */
         ScalarType   mT;
   };

   template <typename TSimType> inline QSTField<TSimType>::QSTField(SmartTruncation pTrunc)
      :SpectralFieldBase(pTrunc), mQ(pTrunc), mS(pTrunc), mT(pTrunc)
   {
   }

   template <typename TSimType> inline const typename QSTField<TSimType>::ScalarType& QSTField<TSimType>::q() const
   {
      return this->mQ;
   }

   template <typename TSimType> inline  typename QSTField<TSimType>::ScalarType& QSTField<TSimType>::rQ()
   {
      return this->mQ;
   }

   template <typename TSimType> inline const  typename QSTField<TSimType>::ScalarType& QSTField<TSimType>::s() const
   {
      return this->mS;
   }

   template <typename TSimType> inline  typename QSTField<TSimType>::ScalarType& QSTField<TSimType>::rS()
   {
      return this->mS;
   }

   template <typename TSimType> inline const  typename QSTField<TSimType>::ScalarType& QSTField<TSimType>::t() const
   {
      return this->mT;
   }

   template <typename TSimType> inline  typename QSTField<TSimType>::ScalarType& QSTField<TSimType>::rT()
   {
      return this->mT;
   }

   template <typename TSimType> void QSTField<TSimType>::computeQSpectra(const typename QSTField<TSimType>::RadialBasisType &radBasis)
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
         this->rQ().rSpectrumL()(l_) += tmpSpectrum.sum();
      }

      // Get the "global" spectra for MPI code
      #ifdef EPMDYNAMO_MPI
         MPI_Allreduce(MPI_IN_PLACE, this->rQ().rSpectrumL().data(), this->q().spectrumL().size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
         MPI_Allreduce(MPI_IN_PLACE, this->rQ().rSpectrumM().data(), this->q().spectrumM().size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      #endif // EPMDYNAMO_MPI
   }

   template <typename TSimType> void QSTField<TSimType>::computeSSpectra(const typename QSTField<TSimType>::RadialBasisType &radBasis)
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
         this->rS().rSpectrumL()(l_) += tmpSpectrum.sum();
      }

      // Get the "global" spectra for MPI code
      #ifdef EPMDYNAMO_MPI
         MPI_Allreduce(MPI_IN_PLACE, this->rS().rSpectrumL().data(), this->s().spectrumL().size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
         MPI_Allreduce(MPI_IN_PLACE, this->rS().rSpectrumM().data(), this->s().spectrumM().size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      #endif // EPMDYNAMO_MPI
   }

   template <typename TSimType> void QSTField<TSimType>::computeTSpectra(const typename QSTField<TSimType>::RadialBasisType &radBasis)
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
         this->rT().rSpectrumL()(l_) += tmpSpectrum.sum();
      }

      // Get the "global" spectra for MPI code
      #ifdef EPMDYNAMO_MPI
         MPI_Allreduce(MPI_IN_PLACE, this->rT().rSpectrumL().data(), this->t().spectrumL().size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
         MPI_Allreduce(MPI_IN_PLACE, this->rT().rSpectrumM().data(), this->t().spectrumM().size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      #endif // EPMDYNAMO_MPI
   }
}

#endif // QSTFIELD_HPP
