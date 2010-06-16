/** \file FDSHOrderScalarBase.hpp
 *  \brief Base of the order ordered FDSH data type implementation
 */

#ifndef FDSHORDERSCALARBASE_HPP
#define FDSHORDERSCALARBASE_HPP

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
    * \brief Base of the order ordered FDSH data type implementation
    */
   class FDSHOrderScalarBase
   {
      public:
         /**
          * @brief Constructor
          *
          * @param pTrunc Truncation information
          */
         FDSHOrderScalarBase(const SmartTruncation pTrunc);

         /**
          * @brief Destructor
          */
         virtual ~FDSHOrderScalarBase() {};

         /**
          * @brief Get the number of harmonic degrees
          *
          * @param m Harmonic order m
          */
         int   nL(const int m) const;

         /**
          * @brief Get the radial size
          *
          * @param m Harmonic order m
          */
         int   nR(const int m) const;

         /**
          * @brief Get the radial sizes
          */
         const ArrayI&  nR() const;

         /**
          * @brief Get the number of harmonic orders
          */
         int   nM() const;

         /**
          * @brief Get the symmetry assumption
          */
         int   mp() const;
         
      protected:
         /**
          * @brief Spectral truncation information
          */
         SmartTruncation   mpTrunc;

         /**
          * @brief Array of nL sizes
          */
         ArrayI mLSizes;

      private:
   };

   inline int FDSHOrderScalarBase::nL(const int m) const
   {
      return this->mLSizes(m);
   }

   inline int FDSHOrderScalarBase::nR(const int m) const
   {
      return this->mpTrunc->local()->fdsh()->nR(m);
   }

   inline const ArrayI&  FDSHOrderScalarBase::nR() const
   {
      return this->mpTrunc->local()->fdsh()->nRArray();
   }

   inline int   FDSHOrderScalarBase::nM() const
   {
      return this->mpTrunc->local()->fdsh()->nM();
   }

   inline int   FDSHOrderScalarBase::mp() const
   {
      return this->mpTrunc->sim()->hoz()->mp();
   }

}

#endif // FDSHORDERSCALARBASE_HPP
