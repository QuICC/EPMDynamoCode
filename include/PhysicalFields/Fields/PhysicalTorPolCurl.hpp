/** \file PhysicalTorPolCurl.hpp
 *  \brief Implementation of Toroidal/Poloidal expanded field with curl
 */

#ifndef PHYSICALTORPOLCURL_HPP
#define PHYSICALTORPOLCURL_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "Simulations/Traits/SimulationTraits.hpp"
#include "GeneralFields/TorPolField.hpp"

namespace EPMDynamo {

   /**
    * \brief Implementation of Toroidal/Poloidal expanded field
    *
    * The used spectral expansion is a Toroidal/Poloidal expansion to take advantage of
    * the divergence free aspect of the field
    *
    * \tparam TSimType Type of the simulation
    * \tparam TBase  Base of the field (used to include imposed field)
    */
   template <typename TSimType, template <typename> class TBase> class PhysicalTorPolCurl: public TBase<TSimType> 
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

   template<typename TSimType, template <typename> class TBase> PhysicalTorPolCurl<TSimType, TBase>::PhysicalTorPolCurl(SmartTruncation pTrunc, typename PhysicalTorPolCurl<TSimType, TBase>::TransformType &transform)
      : TBase<TSimType>(pTrunc, transform)
   {
   }

   template<typename TSimType, template <typename> class TBase> inline void PhysicalTorPolCurl<TSimType, TBase>::transform(const int step)
   {
      if(step == this->mNeedTransform)
      {
         this->updateTotalField();

         this->mrTransform.transformTorPol2RTP(this->rRTP(), this->totalField());

         this->mNeedTransform++;
      }
   }

   template<typename TSimType, template <typename> class TBase> inline void PhysicalTorPolCurl<TSimType, TBase>::curlTransform(const int step)
   {
      if(step == this->mNeedCurlTransform)
      {
         this->updateTotalField();

         this->mrTransform.transformTorPol2Curl(this->rCurl(), this->totalField());

         this->mNeedCurlTransform++;
      }
   }

   template <typename TSimType, template <typename> class TBase> void PhysicalTorPolCurl<TSimType, TBase>::updateSpectra()
   {
      this->rPerturbation().computeTorSpectra(this->mrTransform.radBasis());

      this->rPerturbation().computePolSpectra(this->mrTransform.radBasis());
   }

   template <typename TSimType, template <typename> class TBase> Array PhysicalTorPolCurl<TSimType, TBase>::energy() const
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

   template <typename TSimType, template <typename> class TBase> Matrix PhysicalTorPolCurl<TSimType, TBase>::spectrumL() const
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

   template <typename TSimType, template <typename> class TBase> Matrix PhysicalTorPolCurl<TSimType, TBase>::spectrumM() const
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
