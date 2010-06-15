/** \file SpecTruncation.hpp
 *  \brief Implementation of the remote (core) wide spectral space truncation information
 */

#ifndef SPECTRUNCATION_HPP
#define SPECTRUNCATION_HPP

// Configuration includes
//
#include "Config/SmartPointer.h"

// System includes
//
#include <vector>

// External includes
//

// Project includes
//
#include "General/EPMTypedefs.hpp"

namespace EPMDynamo {

   /**
    * \brief Implementation of the remote (core) wide spectral space truncation information
    */
   class SpecTruncation
   {
      public:
         /**
          * @brief Constructor
          *
          * \param ls Array of harmonic degrees l
          * \param ms Vector of array of harmonic orders m per harmonic degree l
          */
         SpecTruncation(const ArrayI &ls, const std::vector<ArrayI> &ms);

         /**
          * @brief Destructor
          */
         virtual ~SpecTruncation() {};

         /**
          * @brief Get the number of harmonic degrees
          */
         int nL() const;

         /**
          * @brief Get the used harmomic degrees
          */
         const ArrayI& lArray() const;

         /**
          * @brief Get Array of number of harmonic orders
          */
         const ArrayI& nMArray() const;

         /**
          * @brief Get the number of harmonic orders
          *
          * \param l Harmonic degree l
          */
         int nM(const int l) const;

         /**
          * @brief Get the array of used harmonic orders
          *
          * \param l Harmonic degree l
          */
         const ArrayI&  mArray(const int l) const;

         /**
          * @brief Get the used harmonic orders
          *
          * \param idx harmonic order index
          * \param l Harmonic degree l
          */
         int  m(const int idx, const int l) const;
         
      protected:

      private:
         /**
          * @brief The number of harmonic degrees
          */
         int mNl;

         /**
          * @brief The used harmonic degrees
          */
         ArrayI   mLs;

         /**
          * @brief The used harmonic orders
          */
         std::vector<ArrayI>  mMs;

         /**
          * @brief The Array of number of harmonic orders
          */
         ArrayI   mNm;

         /**
          * @brief Initialise the private data
          */
         void init();
   };

   inline int SpecTruncation::nL() const
   {
      return this->mNl;
   }

   inline const ArrayI& SpecTruncation::lArray() const
   {
      return this->mLs;
   }

   inline const ArrayI& SpecTruncation::nMArray() const
   {
      return this->mNm;
   }

   inline int SpecTruncation::nM(const int l) const
   {
      return this->mNm(l);
   }

   inline const ArrayI& SpecTruncation::mArray(const int l) const
   {
      return this->mMs.at(l);
   }

   inline int SpecTruncation::m(const int idx, const int l) const
   {
      return this->mMs.at(l)(idx);
   }

   /// Typedef for a smart SpecTruncation
   typedef EPMSHARED_PTR<SpecTruncation> SmartCoreSpecTrunc;

}

#endif // SPECTRUNCATION_HPP
