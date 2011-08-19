/** \file PhysicalTorPolField.hpp
 *  \brief Implementation of Toroidal/Poloidal expanded field
 */

#ifndef PHYSICALTORPOLFIELD_HPP
#define PHYSICALTORPOLFIELD_HPP

// Configuration includes
//
#include "Config/SimulationConfig.hpp"

// System includes
//

// External includes
//

// Project includes
//
#include "GeneralFields/TorPolField.hpp"

namespace EPMDynamo {

   /**
    * \brief Implementation of Toroidal/Poloidal expanded field
    *
    * The used spectral expansion is a Toroidal/Poloidal expansion to take advantage of
    * the divergence free aspect of the field
    *
    * \tparam TBase  Base of the field (used to include imposed field)
    */
   template <typename TBase> class PhysicalTorPolField: public TBase 
   {
      public:
         /// Typedef from Simulation trait to local scalar type
         typedef SimulationConfig::NumericalScheme::ScalarType    ScalarType;

         /// Typedef from Simulation trait to local transform type
         typedef SimulationConfig::TransformType    TransformType;

         /// Typedef for the spectral field type
         typedef TorPolField  SpectralFieldType;

         /**
          * @brief Constructs the underlying rtp and spectral fields
          *
          * @param pTrunc Truncation information
          * @param transform Reference to the transform object
          */
         PhysicalTorPolField(SmartTruncation pTrunc, TransformType &transform);

         /**
         * @brief Destructor
         */
         virtual ~PhysicalTorPolField() {};

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
          * @brief Generic method to get the symmetric/anti-symmetric energy partition
          */
         virtual Array parityEnergy() const;

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

   template<typename TBase> PhysicalTorPolField<TBase>::PhysicalTorPolField(SmartTruncation pTrunc, typename PhysicalTorPolField<TBase>::TransformType &transform)
      : TBase(pTrunc, transform)
   {
   }

   template<typename TBase> inline void PhysicalTorPolField<TBase>::transform(const int step)
   {
      if(step == this->mNeedTransform)
      {
         this->updateTotalField();

         this->mrTransform.transformTorPol2RTP(this->rRTP(), this->totalField());

         this->mNeedTransform++;
      }
   }

   template <typename TBase> void PhysicalTorPolField<TBase>::updateSpectra()
   {
      this->rPerturbation().computeTorSpectra(this->mrTransform.radBasis());

      this->rPerturbation().computePolSpectra(this->mrTransform.radBasis());
   }

   template <typename TBase> Array PhysicalTorPolField<TBase>::energy() const
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

   template <typename TBase> Array PhysicalTorPolField<TBase>::parityEnergy() const
   {
      // Toroidal symmetric/anti-symmetric energy
      EPMFloat evenTor = this->perturbation().tor().evenEnergy();
      EPMFloat oddTor = this->perturbation().tor().oddEnergy();
      // Poloidal symmetric/anti-symmetric energy
      EPMFloat evenPol = this->perturbation().pol().evenEnergy();
      EPMFloat oddPol = this->perturbation().pol().oddEnergy();

      // Total energy
      Array   energy(6);

      energy(0) = oddTor + evenPol; 
      energy(1) = evenTor + oddPol; 
      energy(2) = oddTor; 
      energy(3) = evenTor; 
      energy(4) = evenPol; 
      energy(5) = oddPol; 

      return energy;
   }

   template <typename TBase> Matrix PhysicalTorPolField<TBase>::spectrumL() const
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

   template <typename TBase> Matrix PhysicalTorPolField<TBase>::spectrumM() const
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

#endif // PHYSICALTORPOLFIELD_HPP
