/** \file PhysicalScalar.hpp
 *  \brief Implementation of physical scalar
 */

#ifndef PHYSICALSCALAR_HPP
#define PHYSICALSCALAR_HPP

// Configuration includes
//
#include "Config/SimulationConfig.hpp"

// System includes
//

// External includes
//

// Project includes
//
#include "General/EPMTypedefs.hpp"

namespace EPMDynamo {

   /**
    * \brief Implementation of physical scalar
    *
    * \tparam TBase  Base of the field (used to include imposed field)
    */
   template <typename TBase> class PhysicalScalar: public TBase
   {
      public:
         /// Typedef from Simulation trait to local scalar type
         typedef SimulationConfig::NumericalScheme::ScalarType  ScalarType;

         /// Typedef from Simulation trait to local transform type
         typedef SimulationConfig::TransformType  TransformType;

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
          * @brief Compute RTP values the scalar
          *
          * @param step Current step in a multistep transform
          */
         virtual void transform(const int step);

         /**
          * @brief Generic method to update the energy spectra
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

         /**
          * @brief Generic method to get the power spectrum per radial polynomial
          */
         virtual Matrix spectrumN() const;
         
      protected:

      private:
   };
   
   template <typename TBase> PhysicalScalar<TBase>::PhysicalScalar(SmartTruncation pTrunc, typename PhysicalScalar<TBase>::TransformType &transform)
      : TBase(pTrunc, transform)
   {
   }

   template <typename TBase> inline void PhysicalScalar<TBase>::transform(const int step)
   {
      if(step == this->mNeedTransform)
      {
         this->updateTotalField();

         this->mrTransform.transformSpec2RTP(this->rRTP(), this->totalField());

         this->mNeedTransform++;
      }
   }

   template <typename TBase> inline void PhysicalScalar<TBase>::updateSpectra()
   {
      this->rPerturbation().computeSpectra(this->mrTransform.radBasis());
   }

   template <typename TBase> inline Array PhysicalScalar<TBase>::energy() const
   {
      Array energy = this->perturbation().energy();

      return energy;
   }

   template <typename TBase> inline Array PhysicalScalar<TBase>::parityEnergy() const
   {
      Array energy(2);
      energy(0) = this->perturbation().evenEnergy();
      energy(1) = this->perturbation().oddEnergy();

      return energy;
   }

   template <typename TBase> inline Matrix PhysicalScalar<TBase>::spectrumN() const
   {
      Matrix   spectrum(this->trunc()->sim()->rad()->nN(),this->trunc()->sim()->hoz()->nL()+1); 

      spectrum.col(0) = Array::LinSpaced(this->trunc()->sim()->rad()->nN(), 0.0, this->trunc()->sim()->rad()->maxN());

      spectrum.rightCols(this->trunc()->sim()->hoz()->nL()) = this->perturbation().spectrumN();

      return spectrum;
   }

   template <typename TBase> inline Matrix PhysicalScalar<TBase>::spectrumL() const
   {
      Matrix   spectrum(this->trunc()->sim()->hoz()->nL(),2); 

      ArrayI lArray = this->trunc()->sim()->hoz()->lArray();

      spectrum.col(0) = lArray.cast<EPMFloat>();

      spectrum.col(1) = this->perturbation().spectrumL();

      return spectrum;
   }

   template <typename TBase> inline Matrix PhysicalScalar<TBase>::spectrumM() const
   {
      Matrix   spectrum(this->trunc()->sim()->hoz()->nM(),2); 
      
      ArrayI mArray = this->trunc()->sim()->hoz()->mArray();

      spectrum.col(0) = mArray.cast<EPMFloat>();

      spectrum.col(1) = this->perturbation().spectrumM();

      return spectrum;
   }

}

#endif // PHYSICALSCALAR_HPP
