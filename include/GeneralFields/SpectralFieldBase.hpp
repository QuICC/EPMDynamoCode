/** \file SpectralFieldBase.hpp
 *  \brief Base for spectral vector field implementations
 */

#ifndef SPECTRALFIELDBASE_HPP
#define SPECTRALFIELDBASE_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "General/EPMTypedefs.hpp"
#include "Domain/Truncation.hpp"

namespace EPMDynamo {

   /**
    * \brief Base for spectral vector field implementations
    */
   class SpectralFieldBase
   {
      public:
         /**
          * @brief Constructor
          *
          * @param pTrunc Truncation information
          */
         SpectralFieldBase(SmartTruncation pTrunc);

         /**
          * @brief Destructor
          */
         virtual ~SpectralFieldBase() {};

         /**
          * @brief Get spectral truncation information pointer
          */
         const SmartTruncation   trunc() const;

         /**
          * @brief Get the number of radial modes
          */
         int nN() const;

         /**
          * @brief Get the number of harmonic degrees
          */
         int nL() const;

         /**
          * @brief Get the number of harmonic orders
          *
          * @param l Harmonic degree
          */
         int nM(const int l) const;
         
      protected:
         /**
          * @brief Pointer to spectral truncation information
          */
         SmartTruncation   mpTrunc;

      private:
   };

   inline const SmartTruncation SpectralFieldBase::trunc() const
   {
      return this->mpTrunc;
   }

   inline int SpectralFieldBase::nN() const
   {
      return this->mpTrunc->sim()->rad()->nN();
   }

   inline int SpectralFieldBase::nL() const
   {
      return this->mpTrunc->local()->spec()->nL();
   }

   inline int SpectralFieldBase::nM(const int l) const
   {
      return this->mpTrunc->local()->spec()->nM(l);
   }

}

#endif // SPECTRALFIELDBASE_HPP
