/** \file PhysicalBase.hpp
 *  \brief Base of the implementation of the physical fields
 */

#ifndef PHYSICALBASE_HPP
#define PHYSICALBASE_HPP

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

namespace EPMDynamo {

   /**
    * \brief Base of the implementation of the physical fields
    */
   class PhysicalBase
   {
      public:
         /// Typedef from Simulation trait to local transform type
         typedef SimulationConfig::TransformType    TransformType;

         /**
         * @brief Construct the general shared information for a physical variable
         *
         * @param pTrunc Truncation information
         * @param transform Spectral transform object reference
         */
         PhysicalBase(SmartTruncation pTrunc, TransformType  &transform);

         /**
         * @brief Destructor
         */
         virtual ~PhysicalBase() {};

         /**
          * @brief Get Spectral truncation information
          */
         const SmartTruncation   trunc() const;

         /**
          * @brief Update the value of the total field (perturbation + imposed)
          */
         virtual void updateTotalField() {};

         /**
          * @brief Generic methode name for perforing the transform to RTP
          *
          * @param step Current step in a multistep transform
          */
         virtual void transform(const int step) = 0;

         /**
          * @brief Generic method to update the energy spectra
          */
         virtual void updateSpectra() = 0;

         /**
          * @brief Generic method to get the energies
          */
         virtual Array energy() const = 0;

         /**
          * @brief Generic method to get the power spectrum per harmonic degree
          */
         virtual Matrix spectrumL() const = 0;

         /**
          * @brief Generic method to get the power spectrum per harmonic order
          */
         virtual Matrix spectrumM() const = 0;
         
      protected:
         /**
          * @brief Reference to Transform object
          */
         TransformType&   mrTransform;

         /**
          * @brief Flag to check if transform has already been computed
          */
         int  mNeedTransform;

      private:

         /**
          * @brief Pointer to spectral truncation information
          */
         SmartTruncation   mpTrunc;
   };

   inline const SmartTruncation PhysicalBase::trunc() const
   {
      return this->mpTrunc;
   }

}

#endif // PHYSICALBASE_HPP
