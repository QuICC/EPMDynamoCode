/** \file LocalTruncation.hpp
 *  \brief Implementation of the local truncation information
 */

#ifndef LOCALTRUNCATION_HPP
#define LOCALTRUNCATION_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "General/EPMTypedefs.hpp"
#include "Domain/Core/RTPDomain.hpp"
#include "Domain/CoreTruncationBase.hpp"

namespace EPMDynamo {

   /**
    * \brief Implementation of the local (current core) truncation information
    */
   class LocalTruncation: public CoreTruncationBase
   {
      public:
         /**
          * @brief Constructor
          *
          * \param rR0 Radial start index
          * \param rNr Number of radial indexes
          * \param rTh0 Array of theta start indexes
          * \param rNth Array of number of theta indexes
          * \param fR0 Array of radial start indexes
          * \param fNr Array of number of radial indexes
          * \param fM Array of harmonic orders
          * \param sL Array of harmonic degrees
          * \param sMs Vector of array of harmonic orders per harmonic degree
          */
         LocalTruncation(const int rR0, const int rNr, const ArrayI &rTh0, const ArrayI &rNth, const ArrayI &fR0, const ArrayI &fNr, const ArrayI &fM, const ArrayI &sL, const std::vector<ArrayI> &sMs);

         /**
          * @brief Destructor
          */
         virtual ~LocalTruncation() {};

         /**
          * @brief Get the RTP truncation
          */
         const SmartCoreRTPDomain rtp() const;

      protected:

      private:
         /**
          * @brief Pointer to the RTP truncation
          */
         SmartCoreRTPDomain mpRTP;
   };

   inline const SmartCoreRTPDomain LocalTruncation::rtp() const
   {
      return this->mpRTP;
   }

   /// Typedef for a smart LocalTruncation
   typedef EPMSHARED_PTR<LocalTruncation> SmartLocalTrunc;

}

#endif // LOCALTRUNCATION_HPP
