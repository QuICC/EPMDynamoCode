/** \file FDSHDegreeScalarBase.hpp
 *  \brief Base of the FDSHDegree ordered data type
 */

#ifndef FDSHDEGREESCALARBASE_HPP
#define FDSHDEGREESCALARBASE_HPP

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
    * \brief Base of the FDSHDegree ordered data type
    *
    * This provides some common methods required for all FDSH fields
    */
   class FDSHDegreeScalarBase
   {
      public:
         /**
         * @brief Stores the pointer to a truncation object
         *
         * @param pTrunc Truncation information
         */
         FDSHDegreeScalarBase(const SmartTruncation pTrunc);

         /**
          * @brief Destructor
          */
         virtual ~FDSHDegreeScalarBase() {};

         /**
          * @brief Get the number of spherical shells
          */
         int   nR() const;

         /**
          * @brief Get the number of harmonic orders
          *
          * @param l Harmonic degree l
          */
         int   nM(const int l) const;

         /**
          * @brief Get the number of harmonic degrees
          */
         int   nL() const;

         /**
          * @brief Get the symmetry assumption
          */
         int   mp() const;

      protected:
         /**
          * @brief Spectral truncation information
          */
         SmartTruncation   mpTrunc;

      private:
   };

   inline int   FDSHDegreeScalarBase::nR() const
   {
      return this->mpTrunc->sim()->rad()->nR();
   }

   inline int   FDSHDegreeScalarBase::nM(const int l) const
   {
      return this->mpTrunc->local()->spec()->nM(l);
   }

   inline int   FDSHDegreeScalarBase::nL() const
   {
      return this->mpTrunc->local()->spec()->nL();
   }

   inline int   FDSHDegreeScalarBase::mp() const
   {
      return this->mpTrunc->sim()->hoz()->mp();
   }

}

#endif // FDSHDEGREESCALARBASE_HPP
