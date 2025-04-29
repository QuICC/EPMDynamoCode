/** \file SpectralSHScalarBase.hpp
 *  \brief Base for a spectral space scalar field implementation
 */

#ifndef SPECTRALSHSCALARBASE_HPP
#define SPECTRALSHSCALARBASE_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "General/EPMTypedefs.hpp"
#include "Domain/Truncation.hpp"
#include "GeneralFields/SpectralFieldBase.hpp"

namespace EPMDynamo {

   /**
    * @brief Basis class for representing a scalar field in full spectral space
    *
    * This datatype is used for storing the full spectral expansion of the scalar field.
    * It will also be the building block of the spectral vector fields.
    */
   class SpectralSHScalarBase: public SpectralFieldBase
   {
      public:
         /**
          * @brief Constructor
          *
          * @param pTrunc truncation information
          * @param needsL0 Is the l=0 mode required in computation?
          */
         SpectralSHScalarBase(const SmartTruncation pTrunc, const bool needsL0);

         /**
          * @brief Destructor
          */
         virtual ~SpectralSHScalarBase() {};

         /**
          * @brief Set the vector of SpectralLShell
          */
         std::vector<SpectralSHLShell>& data();

         /**
          * @brief Get the vector of SpectralSHLShell
          */
         const std::vector<SpectralSHLShell>& data() const;

         /**
          * @brief Get the SpectralSHLShell for a given degree l
          *
          * @param l Requested harmonic degree l
          */
         const SpectralSHLShell&  lshell(const int l) const;
         
         /**
          * @brief Set the SpectralSHLShell for a given degree l
          *
          * @param l Requested harmonic degree l
          */
         SpectralSHLShell&  rLShell(const int l);

         /**
          * @brief Rescale the coefficients by a constant
          *
          * @param scale Scaling factor
          */
         void rescale(const EPMFloat scale);

         /**
          * @brief Lowest harmonic degree l
          */
         int minL() const;
         
      protected:
         /**
          * @brief The lowest harmonic degree (index) stored
          */
         int mMinL;

         /**
          * @brief Storage for the SpectralSHLShell
          */
         std::vector<SpectralSHLShell>    mLShells;

      private:
         /**
          * @brief Initialise storage for the LShells
          */
         void initLShells();
   };

   inline std::vector<SpectralSHLShell>& SpectralSHScalarBase::data()
   {
      return this->mLShells;
   }

   inline const std::vector<SpectralSHLShell>& SpectralSHScalarBase::data() const
   {
      return this->mLShells;
   }

   inline const SpectralSHLShell& SpectralSHScalarBase::lshell(const int l) const
   {
      return this->mLShells.at(l);
   }

   inline SpectralSHLShell& SpectralSHScalarBase::rLShell(const int l)
   {
      return this->mLShells.at(l);
   }

   inline int SpectralSHScalarBase::minL() const
   {
      return this->mMinL;
   }

}

#endif // SPECTRALSHSCALARBASE_HPP
