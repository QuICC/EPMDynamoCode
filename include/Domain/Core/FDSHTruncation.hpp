/** \file FDSHTruncation.hpp
 *  \brief Implementation of the remote (core) FDSH space truncation information
 */

#ifndef FDSHTRUNCATION_HPP
#define FDSHTRUNCATION_HPP

// Configuration includes
//
#include "Config/SmartPointer.h"

// System includes
//

// External includes
//

// Project includes
//
#include "General/EPMTypedefs.hpp"

namespace EPMDynamo {

   /**
    * \brief Implementation of the remote (core) FDSH space truncation information
    */
   class FDSHTruncation
   {
      public:
         /**
          * @brief Constructor
          *
          * \param r0 Array of radial start indexes
          * \param nR Array of number of radial indexes
          * \param mArray Array of harmonic orders
          */
         FDSHTruncation(const ArrayI &r0, const ArrayI &nR, const ArrayI &mArray);

         /**
          * @brief Destructor
          */
         virtual ~FDSHTruncation() {};

         /**
          * @brief Get the number of radial grid points
          *
          * \param m Harmonic order m
          */
         int nR(const int m) const;

         /**
          * @brief Get the number of radial grid points
          */
         const ArrayI& nRArray() const;

         /**
          * @brief Get the radial index of the first radial grid point
          *
          * \param m Harmonic order m
          */
         int r0(const int m) const;

         /**
          * @brief Get the radial index of the first radial grid point
          */
         const ArrayI& r0Array() const;

         /**
          * @brief Get the number of harmonic orders
          */
         int nM() const;

         /**
          * @brief Get the array of harmonic orders
          */
         const ArrayI&   mArray() const;

         /**
          * @brief Get the array of number of harmonic degrees per order
          *
          * \param nL Number of harmonic degrees l
          */
         ArrayI  nLArray(const int nL) const;
         
      protected:

      private:
         /**
          * @brief Index of first radial grid point (from smulation)
          */
         ArrayI mR0;

         /**
          * @brief The number of radial grid points
          */
         ArrayI mNr;

         /**
          * @brief The number of harmonic orders
          */
         int mNm;

         /**
          * @brief The used harmonic orders
          */
         ArrayI   mMArray;
   };

   inline int FDSHTruncation::nR(const int m) const
   {
      return this->mNr(m);
   }

   inline const ArrayI& FDSHTruncation::nRArray() const
   {
      return this->mNr;
   }

   inline int FDSHTruncation::r0(const int m) const
   {
      return this->mR0(m);
   }

   inline const ArrayI& FDSHTruncation::r0Array() const
   {
      return this->mR0;
   }

   inline const ArrayI& FDSHTruncation::mArray() const
   {
      return this->mMArray;
   }

   inline int FDSHTruncation::nM() const
   {
      return this->mNm;
   }

   /// Typedef for a smart FDSHTruncation
   typedef EPMSHARED_PTR<FDSHTruncation> SmartCoreFDSHTrunc;

}

#endif // FDSHTRUNCATION_HPP
