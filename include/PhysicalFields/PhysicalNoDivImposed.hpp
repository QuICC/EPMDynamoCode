/** \file PhysicalNoDivImposed.hpp
 *  \brief Implementation of divergence free field with imposed field
 */

#ifndef PHYSICALNODIVIMPOSED_HPP
#define PHYSICALNODIVIMPOSED_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "Simulations/SimulationTraits.hpp"
#include "PhysicalFields/PhysicalNoDivBase.hpp"

namespace EPMDynamo {

   /**
    * \brief Implementation of divergence free field with imposed field
    *
    * \tparam TSimType Type of the simulation
    */
   template <typename TSimType> class PhysicalNoDivImposed : public PhysicalNoDivBase<TSimType>
   {
      public:
         /// Typedef from Simulation trait to local transform type
         typedef typename SimulationTraits<TSimType>::TransformType    TransformType;

         /// Typedef for the spectral field type
         typedef typename PhysicalNoDivBase<TSimType>::SpectralFieldType  SpectralFieldType;

         /**
          * @brief Constructs the underlying rtp and spectral fields
          *
          * @param pTrunc Truncation information
          * @param transform Reference to the transform object
          */
         PhysicalNoDivImposed(SmartTruncation pTrunc, TransformType &transform);

         /**
         * @brief Destructor
         */
         virtual ~PhysicalNoDivImposed() {};

         /**
          * @brief Get decomposition of the field (total field)
          */
         const SpectralFieldType&  totalField() const;

         /**
          * @brief Get imposed field
          */
         const SpectralFieldType&  imposed() const;

         /**
          * @brief set imposed field 
          */
         SpectralFieldType&  rImposed();

         /**
          * @brief Update the value of the total field (perturbation + imposed)
          */
         void updateTotalField();

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
         /**
          * @brief Spectral decomposition of the field (total field)
          */
         SpectralFieldType    mTotalField;

         /**
          * @brief Spectral decomposition of the imposed field
          */
         SpectralFieldType    mImposedField;

      private:
   };

   template<typename TSimType> PhysicalNoDivImposed<TSimType>::PhysicalNoDivImposed(SmartTruncation pTrunc, typename PhysicalNoDivImposed<TSimType>::TransformType &transform)
      : PhysicalNoDivBase<TSimType>(pTrunc, transform), mTotalField(pTrunc), mImposedField(pTrunc)
   {
   }

   template<typename TSimType> inline const typename PhysicalNoDivImposed<TSimType>::SpectralFieldType& PhysicalNoDivImposed<TSimType>::totalField() const
   {
      return this->mTotalField;
   }

   template<typename TSimType> inline const typename PhysicalNoDivImposed<TSimType>::SpectralFieldType& PhysicalNoDivImposed<TSimType>::imposed() const
   {
      return this->mImposedField;
   }

   template<typename TSimType> inline typename PhysicalNoDivImposed<TSimType>::SpectralFieldType& PhysicalNoDivImposed<TSimType>::rImposed()
   {
      return this->mImposedField;
   }

   template<typename TSimType> inline void PhysicalNoDivImposed<TSimType>::updateTotalField()
   {
      if((this->mNeedTransform == 0) && (this->mNeedCurlTransform == 0))
      {
         // Set toroidal part
         int nL = this->perturbation().tor().nL();
         for(int l = this->perturbation().tor().minL(); l < nL ; ++l)
         {
            this->mTotalField.rTor().rLShell(l) = this->perturbation().tor().lshell(l) + this->imposed().tor().lshell(l);
         }

         // Set poloidal part
         nL = this->perturbation().pol().nL();
         for(int l = this->perturbation().pol().minL(); l < nL ; ++l)
         {
            this->mTotalField.rPol().rLShell(l) = this->perturbation().pol().lshell(l) + this->imposed().pol().lshell(l);
         }
      }
   }

   template<typename TSimType> inline void PhysicalNoDivImposed<TSimType>::transform(const int step)
   {
      if(step == this->mNeedTransform)
      {
         this->updateTotalField();

         this->mrTransform.transformTorPol2RTP(this->rRTP(), this->totalField());

         this->mNeedTransform++;
      }
   }

   template<typename TSimType> inline void PhysicalNoDivImposed<TSimType>::curlTransform(const int step)
   {
      if(step == this->mNeedCurlTransform)
      {
         this->updateTotalField();

         this->mrTransform.transformTorPol2Curl(this->rCurl(), this->totalField());

         this->mNeedCurlTransform++;
      }
   }

   template <typename TSimType> void PhysicalNoDivImposed<TSimType>::updateSpectra()
   {
      this->rPerturbation().computeTorSpectra(this->mrTransform.radBasis());

      this->rPerturbation().computePolSpectra(this->mrTransform.radBasis());
   }

   template <typename TSimType> Array PhysicalNoDivImposed<TSimType>::energy() const
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

   template <typename TSimType> Matrix PhysicalNoDivImposed<TSimType>::spectrumL() const
   {
      // Toroidal Energy spectrum
      Array torE = this->perturbation().tor().spectrumL();
      // Poloidal Energy spectrum
      Array polE = this->perturbation().pol().spectrumL();

      Matrix   spectrum(this->trunc()->sim()->hoz()->nL(), 4); 

      ArrayI lArray = this->trunc()->sim()->hoz()->lArray();
      spectrum.col(0) = lArray.cast<DynamoFloat>();
      
      spectrum.col(1) = torE + polE;
      spectrum.col(2) = torE;
      spectrum.col(3) = polE;

      return spectrum;
   }

   template <typename TSimType> Matrix PhysicalNoDivImposed<TSimType>::spectrumM() const
   {
      // Toroidal Energy spectrum
      Array torE = this->perturbation().tor().spectrumM();
      // Poloidal Energy spectrum
      Array polE = this->perturbation().pol().spectrumM();

      Matrix   spectrum(this->trunc()->sim()->hoz()->nM(), 4); 

      ArrayI mArray = this->trunc()->sim()->hoz()->mArray();
      spectrum.col(0) = mArray.cast<DynamoFloat>();
      
      spectrum.col(1) = torE + polE;
      spectrum.col(2) = torE;
      spectrum.col(3) = polE;

      return spectrum;
   }

}

#endif // PHYSICALNODIVIMPOSED_HPP
