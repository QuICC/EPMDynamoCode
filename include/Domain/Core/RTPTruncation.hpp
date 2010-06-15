/** \file RTPTruncation.hpp
 *  \brief Implementation of the remote (core) wide RTP space truncation information
 */

#ifndef RTPTRUNCATION_HPP
#define RTPTRUNCATION_HPP

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
    * \brief Implementation of the remote (core) wide RTP space truncation information
    */
   class RTPTruncation
   {
      public:
         /**
          * @brief Constructor 
          *
          * \param r0 Radial start index
          * \param nR Number of radial indexes
          * \param th0 Array of theta start indexes
          * \param nTh Array of number of theta values
          */
         RTPTruncation(const int r0, const int nR, const ArrayI &th0, const ArrayI &nTh);

         /**
          * @brief Destructor
          */
         virtual ~RTPTruncation() {};

         /**
          * @brief Get the number of radial grid points
          */
         int nR() const;

         /**
          * @brief Get the radial index of the first radial grid point
          */
         int r0() const;

         /**
          * @brief Get the number of theta grid points for given radial grid
          *
          * @param n Radial grid index
          */
         int nTh(const int n) const;

         /**
          * @brief Get Array of theta grid sizes
          */
         const ArrayI& nThArray() const;

         /**
          * @brief Get the maximum number of theta grid points
          */
         int maxNth() const;

         /**
          * @brief first index of the theta grid points
          *
          * @param n Radial grid index
          */
         int th0(const int n) const;

         /**
          * @brief Array of first index of the theta grid points
          */
         const ArrayI& th0Array() const;
         
      protected:

      private:
         /**
          * @brief Index of first radial grid point (from smulation)
          */
         int mR0;

         /**
          * @brief The number of radial grid points
          */
         int mNr;

         /**
          * @brief Array of first indexes of Theta grid
          */
         ArrayI mTh0;

         /**
          * @brief Array of Theta sizes
          */
         ArrayI mNth;
   };

   inline int RTPTruncation::nR() const
   {
      return this->mNr;
   }

   inline int RTPTruncation::nTh(const int n) const
   {
      return this->mNth(n);
   }

   inline const ArrayI& RTPTruncation::nThArray() const
   {
      return this->mNth;
   }

   inline int RTPTruncation::maxNth() const
   {
      return this->mNth.maxCoeff();
   }

   inline int RTPTruncation::r0() const
   {
      return this->mR0;
   }

   inline int RTPTruncation::th0(const int n) const
   {
      return this->mTh0(n);
   }

   inline const ArrayI& RTPTruncation::th0Array() const
   {
      return this->mTh0;
   }

   /// Typedef for a smart RTPTruncation
   typedef EPMSHARED_PTR<RTPTruncation> SmartCoreRTPTrunc;

}

#endif // RTPTRUNCATION_HPP
