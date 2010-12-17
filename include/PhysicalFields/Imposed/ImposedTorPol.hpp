/** \file ImposedTorPol.hpp
 *  \brief Implementation of Toroidal/Poloidal expanded imposed field
 */

#ifndef IMPOSEDTORPOL_HPP
#define IMPOSEDTORPOL_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "Simulations/Traits/SimulationTraits.hpp"
#include "PhysicalFields/Fields/PhysicalTorPolBase.hpp"

namespace EPMDynamo {

   /**
    * \brief Implementation of Toroidal/Poloidal expanded imposed field
    *
    * \tparam TSimType Type of the simulation
    */
   template <typename TSimType> class ImposedTorPol : public PhysicalTorPolBase<TSimType>
   {
      public:
         /// Typedef from Simulation trait to local transform type
         typedef typename SimulationTraits<TSimType>::TransformType    TransformType;

         /// Typedef for the spectral field type
         typedef typename PhysicalTorPolBase<TSimType>::SpectralFieldType  SpectralFieldType;

         /**
          * @brief Constructs the underlying rtp and spectral fields
          *
          * @param pTrunc Truncation information
          * @param transform Reference to the transform object
          */
         ImposedTorPol(SmartTruncation pTrunc, TransformType &transform);

         /**
         * @brief Destructor
         */
         virtual ~ImposedTorPol() {};

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
         virtual void updateTotalField();
         
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

   template<typename TSimType> ImposedTorPol<TSimType>::ImposedTorPol(SmartTruncation pTrunc, typename ImposedTorPol<TSimType>::TransformType &transform)
      : PhysicalTorPolBase<TSimType>(pTrunc, transform), mTotalField(pTrunc), mImposedField(pTrunc)
   {
   }

   template<typename TSimType> inline const typename ImposedTorPol<TSimType>::SpectralFieldType& ImposedTorPol<TSimType>::totalField() const
   {
      return this->mTotalField;
   }

   template<typename TSimType> inline const typename ImposedTorPol<TSimType>::SpectralFieldType& ImposedTorPol<TSimType>::imposed() const
   {
      return this->mImposedField;
   }

   template<typename TSimType> inline typename ImposedTorPol<TSimType>::SpectralFieldType& ImposedTorPol<TSimType>::rImposed()
   {
      return this->mImposedField;
   }

   template<typename TSimType> inline void ImposedTorPol<TSimType>::updateTotalField()
   {
      if(this->needAnyTransform())
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

}

#endif // IMPOSEDTORPOL_HPP
