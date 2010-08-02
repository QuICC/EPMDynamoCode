/** \file PhysicalScalar.hpp
 *  \brief Implementation of physical scalar
 */

#ifndef PHYSICALSCALAR_HPP
#define PHYSICALSCALAR_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "Simulations/Traits/SimulationTraits.hpp"
#include "General/EPMTypedefs.hpp"

namespace EPMDynamo {

   /**
    * \brief Implementation of physical scalar
    *
    * \tparam TSimType Type of the simulation
    * \tparam TBase  Base of the field (used to include imposed field)
    */
   template <typename TSimType, template <typename> class TBase> class PhysicalScalar: public TBase<TSimType>
   {
      public:
         /// Typedef from Simulation trait to local scalar type
         typedef typename TSimType::ScalarType  ScalarType;

         /// Typedef from Simulation trait to local transform type
         typedef typename SimulationTraits<TSimType>::TransformType  TransformType;

         /// Typedef for the spectral field type
         typedef ScalarType SpectralFieldType;

         /**
         * @brief Constructor
         *
         * @param pTrunc Truncation information
         * @param transform Reference to transform object
         */
         PhysicalScalar(SmartTruncation pTrunc, TransformType &transform);

         /**
         * @brief Destructor
         */
         virtual ~PhysicalScalar() {};

         /**
          * @brief Get spectral expansion scalar (perturbation part)
          */
         const ScalarType&  perturbation() const;

         /**
          * @brief Get spectral expansion scalar (Total scalar field)
          */
         const ScalarType&  totalField() const;

         /**
          * @brief Set spectral expansion scalar
          */
         ScalarType&  rPerturbation();

         /**
          * @brief Compute RTP values the scalar
          *
          * @param step Current step in a multistep transform
          */
         virtual void transform(const int step);

         /**
          * @brief Compute RTP values of the gradient of the scalar
          *
          * @param step Current step in a multistep transform
          */
         virtual void gradTransform(const int step);

         /**
          * @brief Generic method to update the energy spectra
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

         /**
          * @brief Spectral expansion scalar
          */
         ScalarType   mPerturbation;

         /**
          * @brief Flag to check if grad transform has already been computed
          */
         int  mNeedGradTransform;
   };
   
   template <typename TSimType, template <typename> class TBase> PhysicalScalar<TSimType, TBase>::PhysicalScalar(SmartTruncation pTrunc, typename PhysicalScalar<TSimType, TBase>::TransformType &transform)
      : TBase<TSimType>(pTrunc, transform), mPerturbation(pTrunc, true), mNeedGradTransform(0)
   {
   }

   template <typename TSimType, template <typename> class TBase> inline const typename PhysicalScalar<TSimType, TBase>::ScalarType& PhysicalScalar<TSimType, TBase>::perturbation() const
   {
      return this->mPerturbation;
   }

   template <typename TSimType, template <typename> class TBase> inline const typename PhysicalScalar<TSimType, TBase>::ScalarType& PhysicalScalar<TSimType, TBase>::totalField() const
   {
      return this->mPerturbation;
   }

   template <typename TSimType, template <typename> class TBase> inline typename PhysicalScalar<TSimType, TBase>::ScalarType& PhysicalScalar<TSimType, TBase>::rPerturbation()
   {
      this->mNeedTransform = 0;

      this->mNeedGradTransform = 0;

      return this->mPerturbation;
   }

   template <typename TSimType, template <typename> class TBase> inline void PhysicalScalar<TSimType, TBase>::transform(const int step)
   {
      if(step == this->mNeedTransform)
      {
         this->updateTotalField();

         this->mrTransform.transformSpec2RTP(this->rRTP(), this->totalField());

         this->mNeedTransform++;
      }
   }

   template <typename TSimType, template <typename> class TBase> inline void PhysicalScalar<TSimType, TBase>::gradTransform(const int step)
   {
      if(step == this->mNeedGradTransform)
      {
         this->updateTotalField();

         this->mrTransform.transformSpec2Grad(this->rGrad(), this->totalField());

         this->mNeedGradTransform++;
      }
   }

   template <typename TSimType, template <typename> class TBase> inline void PhysicalScalar<TSimType, TBase>::updateSpectra()
   {
      this->rPerturbation().computeSpectra(this->mrTransform.radBasis());
   }

   template <typename TSimType, template <typename> class TBase> inline Array PhysicalScalar<TSimType, TBase>::energy() const
   {
      Array energy = this->perturbation().energy();

      return energy;
   }

   template <typename TSimType, template <typename> class TBase> inline Matrix PhysicalScalar<TSimType, TBase>::spectrumL() const
   {
      Matrix   spectrum(this->trunc()->sim()->hoz()->nL(),2); 

      ArrayI lArray = this->trunc()->sim()->hoz()->lArray();

      spectrum.col(0) = lArray.cast<EPMFloat>();

      spectrum.col(1) = this->perturbation().spectrumL();

      return spectrum;
   }

   template <typename TSimType, template <typename> class TBase> inline Matrix PhysicalScalar<TSimType, TBase>::spectrumM() const
   {
      Matrix   spectrum(this->trunc()->sim()->hoz()->nM(),2); 
      
      ArrayI mArray = this->trunc()->sim()->hoz()->mArray();

      spectrum.col(0) = mArray.cast<EPMFloat>();

      spectrum.col(1) = this->perturbation().spectrumM();

      return spectrum;
   }

}

#endif // PHYSICALSCALAR_HPP
