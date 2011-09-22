/** \file PhysicalTorPolCurl.hpp
 *  \brief Implementation of Toroidal/Poloidal expanded field with curl
 */

#ifndef PHYSICALTORPOLCURL_HPP
#define PHYSICALTORPOLCURL_HPP

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
   template <typename TBase> class PhysicalTorPolCurl: public TBase 
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
         PhysicalTorPolCurl(SmartTruncation pTrunc, TransformType &transform);

         /**
         * @brief Destructor
         */
         virtual ~PhysicalTorPolCurl() {};

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
          * @brief Generic method to get the symmetric/anti-symmetric partition of energy
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

         /**
          * @brief Generic method to get the power spectrum per radial poynomial
          */
         virtual Matrix spectrumN() const;
         
      protected:

      private:
   };

   template<typename TBase> PhysicalTorPolCurl<TBase>::PhysicalTorPolCurl(SmartTruncation pTrunc, typename PhysicalTorPolCurl<TBase>::TransformType &transform)
      : TBase(pTrunc, transform)
   {
   }

   template<typename TBase> inline void PhysicalTorPolCurl<TBase>::transform(const int step)
   {
      if(step == this->mNeedTransform)
      {
         this->updateTotalField();

         this->mrTransform.transformTorPol2RTP(this->rRTP(), this->totalField());

         this->mNeedTransform++;
      }
   }

   template<typename TBase> inline void PhysicalTorPolCurl<TBase>::curlTransform(const int step)
   {
      if(step == this->mNeedCurlTransform)
      {
         this->updateTotalField();

         this->mrTransform.transformTorPol2Curl(this->rCurl(), this->totalField());

         this->mNeedCurlTransform++;
      }
   }

   template <typename TBase> void PhysicalTorPolCurl<TBase>::updateSpectra()
   {
      this->rPerturbation().computeTorSpectra(this->mrTransform.radBasis());

      this->rPerturbation().computePolSpectra(this->mrTransform.radBasis());
   }

   template <typename TBase> Array PhysicalTorPolCurl<TBase>::energy() const
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

   template <typename TBase> Array PhysicalTorPolCurl<TBase>::parityEnergy() const
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

   template <typename TBase> inline Matrix PhysicalTorPolCurl<TBase>::spectrumN() const
   {
      Matrix   spectrum(this->trunc()->sim()->rad()->nN(),this->trunc()->sim()->hoz()->nL()+1); 

      spectrum.col(0) = Array::LinSpaced(this->trunc()->sim()->rad()->nN(), 0.0, this->trunc()->sim()->rad()->maxN());

      spectrum.rightCols(this->trunc()->sim()->hoz()->nL()) = this->perturbation().tor().spectrumN()+this->perturbation().pol().spectrumN();

      return spectrum;
   }

   template <typename TBase> Matrix PhysicalTorPolCurl<TBase>::spectrumL() const
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

   template <typename TBase> Matrix PhysicalTorPolCurl<TBase>::spectrumM() const
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

#endif // PHYSICALTORPOLCURL_HPP
