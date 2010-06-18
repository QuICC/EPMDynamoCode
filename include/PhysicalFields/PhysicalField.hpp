/** \file PhysicalField.hpp
 *  \brief Implementation of a general physical field (QST)
 */

#ifndef PHYSICALFIELD_HPP
#define PHYSICALFIELD_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "Simulations/Traits/SimulationTraits.hpp"
#include "PhysicalFields/PhysicalFieldBase.hpp"
#include "Domain/Truncation.hpp"
#include "GeneralFields/QSTField.hpp"

namespace EPMDynamo {

   /**
    * @brief This class template implements a general physical field
    *
    * The used spectral expansion is an QST expansion as it is a completely general
    * expansion
    *
    * \tparam TSimType Type of the simulation
    */
   template <typename TSimType> class PhysicalField : public PhysicalFieldBase<TSimType>
   {
      public:
         /// Typedef from Simulation trait to local scalar type
         typedef typename TSimType::ScalarType    ScalarType;

         /// Typedef from Simulation trait to local transform type
         typedef typename SimulationTraits<TSimType>::TransformType    TransformType;

         /// Typedef for the spectral field type
         typedef QSTField<TSimType>  SpectralFieldType;

         /**
         * @brief Constructs underlying rtp and spectral fields
         *
         * @param pTrunc Truncation information
         * @param transform Transform object
         */
         PhysicalField(SmartTruncation pTrunc, TransformType &transform);

         /**
         * @brief Simple empty destructor
         */
         virtual ~PhysicalField() {};

         /**
          * @brief Get QST decomposition of the perturbation part
          */
         const QSTField<TSimType>&  perturbation() const;

         /**
          * @brief Get QST decomposition of the total field
          */
         const QSTField<TSimType>&  totalField() const;

         /**
          * @brief Set QST decomposition of the perturbation part
          */
         QSTField<TSimType>&  rPerturbation();

         /**
          * @brief Compute RTP values of the field
          *
          * @param step Current step in a multistep transform
          */
         virtual void transform(const int step);

         /**
          * @brief Compute RTP values of the curl of the field
          *
          * @param step Current step in a multistep transform
          */
         virtual void curlTransform(const int step);

         /**
          * @brief Generic method to get the energies
          */
         virtual Array energy() const;

         /**
          * @brief Generic method to get the power spectrum per harmonic order
          */
         virtual Matrix spectrumM() const;

         /**
          * @brief Generic method to get the power spectrum per harmonic degree
          */
         virtual Matrix spectrumL() const;
         
      protected:

         /**
          * @brief Spectral QST decomposition of the field
          */
         QSTField<TSimType>    mPerturbation;

      private:
         /**
          * @brief Flag to check if curl transform has already been
          *          computed
          */
         int  mNeedCurlTransform;
   };

   template<typename TSimType> PhysicalField<TSimType>::PhysicalField(SmartTruncation pTrunc, typename PhysicalField<TSimType>::TransformType &transform)
      : PhysicalFieldBase<TSimType>(pTrunc, transform), mPerturbation(pTrunc), mNeedCurlTransform(0)
   {
   }

   template<typename TSimType> inline const QSTField<TSimType>& PhysicalField<TSimType>::perturbation() const
   {
      return this->mPerturbation;
   }

   template<typename TSimType> inline const QSTField<TSimType>& PhysicalField<TSimType>::totalField() const
   {
      return this->mPerturbation;
   }

   template<typename TSimType> inline QSTField<TSimType>& PhysicalField<TSimType>::rPerturbation()
   {
      this->mNeedTransform = 0;

      this->mNeedCurlTransform = 0;
      
      return mPerturbation;
   }

   template<typename TSimType> inline void PhysicalField<TSimType>::transform(const int step)
   {
      if(step == this->mNeedTransform)
      {
         this->mrTransform.transformQST2RTP(this->rRTP(), this->totalField());

         this->mNeedTransform++;
      }
   }

   template<typename TSimType> inline void PhysicalField<TSimType>::curlTransform(const int step)
   {
      if(step == this->mNeedCurlTransform)
      {
         this->mrTransform.transformQST2Curl(this->rCurl(), this->totalField());

         this->mNeedCurlTransform++;
      }
   }

   template <typename TSimType> inline void PhysicalField<TSimType>::updateSpectra()
   {
      this->rPerturbation().computeQSpectra(this->mrTransform.radBasis());

      this->rPerturbation().computeSSpectra(this->mrTransform.radBasis());

      this->rPerturbation().computeTSpectra(this->mrTransform.radBasis());
   }

   template <typename TSimType> inline Array PhysicalField<TSimType>::energy() const
   {
      // Q Energy
      Array qE = this->perturbation().q().Energy();
      // S Energy
      Array sE = this->perturbation().s().Energy();
      // T Energy
      Array tE = this->perturbation().t().Energy();

      // Total energy
      Array   energy(4);

      energy(0) = qE(0) + sE(0) + tE(0); 
      energy(1) = qE(0); 
      energy(2) = sE(0); 
      energy(3) = tE(0); 

      return energy;
   }

   template <typename TSimType> inline Matrix PhysicalField<TSimType>::spectrumL() const
   {
      // QST Q Energy spectrum
      Array qE = this->perturbation().q().SpectrumL();
      // QST S Energy spectrum
      Array sE = this->perturbation().s().SpectrumL();
      // QST T Energy spectrum
      Array tE = this->perturbation().t().SpectrumL();

      Matrix   spectrum(this->trunc()->local()->spec()->nL(), 5); 

      ArrayI lArray = this->trunc()->local()->spec()->lArray();
      spectrum.col(0) = lArray.cast<EPMFloat>();
      
      spectrum.col(1) = qE + sE + tE;
      spectrum.col(2) = qE;
      spectrum.col(3) = sE;
      spectrum.col(3) = tE;

      return spectrum;
   }

   template <typename TSimType> inline Matrix PhysicalField<TSimType>::spectrumM() const
   {
      // QST Q Energy spectrum
      Array qE = this->perturbation().q().SpectrumM();
      // QST S Energy spectrum
      Array sE = this->perturbation().s().SpectrumM();
      // QST T Energy spectrum
      Array tE = this->perturbation().t().SpectrumM();

      Matrix   spectrum(this->trunc()->sim()->hoz()->nM(), 5); 

      ArrayI mArray = this->trunc()->sim()->hoz()->mArray();
      spectrum.col(0) = mArray.cast<EPMFloat>();
      
      spectrum.col(1) = qE + sE + tE;
      spectrum.col(2) = qE;
      spectrum.col(3) = sE;
      spectrum.col(3) = tE;

      return spectrum;
   }

}

#endif // PHYSICALFIELD_HPP
