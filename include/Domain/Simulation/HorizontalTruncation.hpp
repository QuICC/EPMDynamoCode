/** \file HorizontalTruncation.hpp
 *  \brief Implementation of the simulation wide horizontal truncation information
 */

#ifndef HORIZONTALTRUNCATION_HPP
#define HORIZONTALTRUNCATION_HPP

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
     * \brief Implementation of the simulation wide horizontal truncation information
     */
   class HorizontalTruncation
   {
      public:
         /**
          * @brief Constructor
          *
          * @param maxL  Highest harmonic degree
          * @param maxM  Highest harmonic order
          * @param mp    Longitudinal symmetry  assumption
          */
         HorizontalTruncation(const int maxL, const int maxM, const int mp);

         /**
          * @brief Destructor
          */
         virtual ~HorizontalTruncation() {};

         /**
          * @brief Get Highest harmonic degree
          */
         int maxL() const;

         /**
          * @brief Get Highest harmonic order
          */
         int maxM() const;

         /**
          * @brief Get the longitudinal symmetry assumption factor
          */
         int mp() const;

         /**
          * @brief Get the total number of harmonics
          */
         int nH() const;

         /**
          * @brief Get the total number of harmonic degrees
          */
         int nL() const;

         /**
          * @brief Get the total number of harmonic orders
          */
         int nM() const;

         /**
          * @brief Get the number of harmonic orders for a given order
          *
          * \param l Harmonic degree l
          */
         int nM(const int l) const;

         /**
          * @brief Get the size of the theta grid
          */
         int nTh() const;

         /**
          * @brief Get the size of the phi grid
          */
         int nPh() const;

         /**
          * @brief Get the array of harmonic orders
          */
         ArrayI   mArray() const;

         /**
          * @brief Get the array of harmonic degrees
          */
         ArrayI   lArray() const;

         /**
          * @brief Get the sizes of the degrees
          */
         ArrayI   nLArray() const;
         
      protected:

      private:
         /**
          * @brief Highest harmonic order L
          */
         int mMaxL;

         /**
          * @brief Highest harmonic degree m
          */
         int mMaxM;

         /**
          * @brief Longitudinal symmetry Mp
          */
         int mMp;

         /**
          * @brief Effective number of harmonic orders L
          */
         int mNl;

         /**
          * @brief Effective number of harmonic degrees m
          */
         int mNm;

         /**
          * @brief Total number of harmonics
          */
         int mNh;

         /**
          * @brief Number of longitudinal (phi) grid points
          */
         int mNph;

         /**
          * @brief Number of latitudinal (theta) grid points
          */
         int mNth;

         /**
          * @brief Initialise total number of harmonics
          */
         void initH();

         /**
          * @brief Theta dealiasing factor
          */
         static const EPMFloat msThetaDealiasing;

         /**
          * @brief Phi dealiasing factor
          */
         static const EPMFloat msPhiDealiasing;

         /**
          * @brief Correct grid sizes if grid type requires it (ex. CSCS)
          */
         void correctGridSizes();
   };

   inline int HorizontalTruncation::maxL() const
   {
      return this->mMaxL;
   }

   inline int HorizontalTruncation::maxM() const
   {
      return this->mMaxM;
   }

   inline int HorizontalTruncation::nH() const
   {
      return this->mNh;
   }

   inline int HorizontalTruncation::nM() const
   {
      return this->mNm;
   }

   inline int HorizontalTruncation::nM(const int l) const
   {
      if(l < this->maxM())
      {
         return l/this->mp() + 1;
      } else 
      {
         return this->nM();
      }
   }

   inline int HorizontalTruncation::nL() const
   {
      return this->mNl;
   }

   inline int HorizontalTruncation::mp() const
   {
      return this->mMp;
   }

   inline int HorizontalTruncation::nTh() const
   {
      return this->mNth;
   }

   inline int HorizontalTruncation::nPh() const
   {
      return this->mNph;
   }

   /// Typedef for a Smart pointer to a HorizontalTruncation
   typedef EPMSHARED_PTR<HorizontalTruncation>  SmartSimHozTrunc;

}

#endif // HORIZONTALTRUNCATION_HPP
