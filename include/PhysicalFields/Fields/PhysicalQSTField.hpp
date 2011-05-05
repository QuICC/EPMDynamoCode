/** \file PhysicalQSTField.hpp
 *  \brief Implementation of QST expanded field
 */

#ifndef PHYSICALQSTFIELD_HPP
#define PHYSICALQSTFIELD_HPP

// Configuration includes
//
#include "Config/SimulationConfig.hpp"

// System includes
//

// External includes
//

// Project includes
//
#include "Domain/Truncation.hpp"
#include "GeneralFields/QSTField.hpp"

namespace EPMDynamo {

   /**
    * \brief Implementation of QST expanded field
    *
    * The used spectral expansion is an QST expansion as it is a completely general
    * expansion
    *
    * \tparam TBase  Base of the field (used to include imposed field)
    */
   template <typename TBase> class PhysicalQSTField : public TBase
   {
      public:
         /// Typedef from Simulation trait to local scalar type
         typedef SimulationConfig::NumericalScheme::ScalarType    ScalarType;

         /// Typedef from Simulation trait to local transform type
         typedef SimulationConfig::TransformType    TransformType;

         /// Typedef for the spectral field type
         typedef QSTField  SpectralFieldType;

         /**
         * @brief Constructs underlying rtp and spectral fields
         *
         * @param pTrunc Truncation information
         * @param transform Transform object
         */
         PhysicalQSTField(SmartTruncation pTrunc, TransformType &transform);

         /**
         * @brief Simple empty destructor
         */
         virtual ~PhysicalQSTField() {};

         /**
          * @brief Compute RTP values of the field
          *
          * @param step Current step in a multistep transform
          */
         virtual void transform(const int step);

         /**
          * @brief Generic method to get the energies
          */
         virtual void updateSpectra();

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

      private:
   };

   template<typename TBase> PhysicalQSTField<TBase>::PhysicalQSTField(SmartTruncation pTrunc, typename PhysicalQSTField<TBase>::TransformType &transform)
      : TBase(pTrunc, transform)
   {
   }

   template<typename TBase> inline void PhysicalQSTField<TBase>::transform(const int step)
   {
      if(step == this->mNeedTransform)
      {
         this->updateTotalField();

         this->mrTransform.transformQST2RTP(this->rRTP(), this->totalField());

         this->mNeedTransform++;
      }
   }

   template <typename TBase> inline void PhysicalQSTField<TBase>::updateSpectra()
   {
      this->rPerturbation().computeQSpectra(this->mrTransform.radBasis());

      this->rPerturbation().computeSSpectra(this->mrTransform.radBasis());

      this->rPerturbation().computeTSpectra(this->mrTransform.radBasis());
   }

   template <typename TBase> inline Array PhysicalQSTField<TBase>::energy() const
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

   template <typename TBase> inline Matrix PhysicalQSTField<TBase>::spectrumL() const
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

   template <typename TBase> inline Matrix PhysicalQSTField<TBase>::spectrumM() const
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

#endif // PHYSICALQSTFIELD_HPP
