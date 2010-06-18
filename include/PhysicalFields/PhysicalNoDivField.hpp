/** \file PhysicalNoDivField.hpp
 *  \brief Implementation of divergence free field
 */

#ifndef PHYSICALNODIVFIELD_HPP
#define PHYSICALNODIVFIELD_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "Simulations/Traits/SimulationTraits.hpp"
#include "PhysicalFields/PhysicalNoDivBase.hpp"
#include "GeneralFields/TorPolField.hpp"

namespace EPMDynamo {

   /**
    * \brief Implementation of divergence free field
    *
    * The used spectral expansion is a Toroidal/Poloidal expansion to take advantage of
    * the divergence free aspect of the field
    *
    * \tparam TSimType Type of the simulation
    */
   template <typename TSimType> class PhysicalNoDivField: public PhysicalNoDivBase<TSimType> 
   {
      public:
         /// Typedef from Simulation trait to local scalar type
         typedef typename TSimType::ScalarType    ScalarType;

         /// Typedef from Simulation trait to local transform type
         typedef typename SimulationTraits<TSimType>::TransformType    TransformType;

         /// Typedef for the spectral field type
         typedef TorPolField<TSimType>  SpectralFieldType;

         /**
          * @brief Constructs the underlying rtp and spectral fields
          *
          * @param pTrunc Truncation information
          * @param transform Reference to the transform object
          */
         PhysicalNoDivField(SmartTruncation pTrunc, TransformType &transform);

         /**
         * @brief Destructor
         */
         virtual ~PhysicalNoDivField() {};

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

   template<typename TSimType> PhysicalNoDivField<TSimType>::PhysicalNoDivField(SmartTruncation pTrunc, typename PhysicalNoDivField<TSimType>::TransformType &transform)
      : PhysicalNoDivBase<TSimType>(pTrunc, transform)
   {
   }

   template<typename TSimType> inline void PhysicalNoDivField<TSimType>::transform(const int step)
   {
      if(step == this->mNeedTransform)
      {
         this->mrTransform.transformTorPol2RTP(this->rRTP(), this->totalField());

         this->mNeedTransform++;
      }
   }

   template<typename TSimType> inline void PhysicalNoDivField<TSimType>::curlTransform(const int step)
   {
      if(step == this->mNeedCurlTransform)
      {
         this->mrTransform.transformTorPol2Curl(this->rCurl(), this->totalField());

         this->mNeedCurlTransform++;
      }
   }

   template <typename TSimType> void PhysicalNoDivField<TSimType>::updateSpectra()
   {
      this->rPerturbation().computeTorSpectra(this->mrTransform.radBasis());

      this->rPerturbation().computePolSpectra(this->mrTransform.radBasis());
   }

   template <typename TSimType> Array PhysicalNoDivField<TSimType>::energy() const
   {
      // Toroidal Energy
      Array torE = this->perturbation().tor().energy();
      // Poloidal Energy
      Array polE = this->perturbation().pol().energy();

      // Total energy
      Array   energy(3);

      energy(0) = torE(0) + polE(0); 
      energy(1) = torE(0); 
      energy(2) = polE(0); 

      return energy;
   }

   template <typename TSimType> Matrix PhysicalNoDivField<TSimType>::spectrumL() const
   {
      // Toroidal Energy spectrum
      Array torE = this->perturbation().tor().spectrumL();
      // Poloidal Energy spectrum
      Array polE = this->perturbation().pol().spectrumL();

      Matrix   spectrum(this->trunc()->sim()->hoz()->nL(), 4); 

      ArrayI lArray = this->trunc()->sim()->hoz()->lArray();
      spectrum.col(0) = lArray.cast<EPMFloat>();
      
      spectrum.col(1) = torE + polE;
      spectrum.col(2) = torE;
      spectrum.col(3) = polE;

      return spectrum;
   }

   template <typename TSimType> Matrix PhysicalNoDivField<TSimType>::spectrumM() const
   {
      // Toroidal Energy spectrum
      Array torE = this->perturbation().tor().spectrumM();
      // Poloidal Energy spectrum
      Array polE = this->perturbation().pol().spectrumM();

      Matrix   spectrum(this->trunc()->sim()->hoz()->nM(), 4); 

      ArrayI mArray = this->trunc()->sim()->hoz()->mArray();
      spectrum.col(0) = mArray.cast<EPMFloat>();
      
      spectrum.col(1) = torE + polE;
      spectrum.col(2) = torE;
      spectrum.col(3) = polE;

      return spectrum;
   }

}

#endif // PHYSICALNODIVFIELD_HPP
