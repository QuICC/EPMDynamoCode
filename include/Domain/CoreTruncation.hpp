/** \file CoreTruncation.hpp
 *  \brief Implementation of the core truncation information
 */

#ifndef CORETRUNCATION_HPP
#define CORETRUNCATION_HPP

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
#include "Domain/Core/RTPTruncation.hpp"
#include "Domain/CoreTruncationBase.hpp"

namespace EPMDynamo {

   /**
    * \brief Implementation of the core (not current one) truncation information
    */
   class CoreTruncation : public CoreTruncationBase
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
         CoreTruncation(const int rR0, const int rNr, const ArrayI &rTh0, const ArrayI &rNth, const ArrayI &fR0, const ArrayI &fNr, const ArrayI &fM, const ArrayI &sL, const std::vector<ArrayI> &sMs);

         /**
          * @brief Destructor
          */
         virtual ~CoreTruncation() {};

         /**
          * @brief Get the RTP truncation
          */
         const SmartCoreRTPTrunc rtp() const;
         
      protected:

      private:
         /**
          * @brief Pointer to the RTP truncation
          */
         SmartCoreRTPTrunc mpRTP;
   };

   inline const SmartCoreRTPTrunc CoreTruncation::rtp() const
   {
      return this->mpRTP;
   }

   /// Typedef for a smart CoreTruncation
   typedef EPMSHARED_PTR<CoreTruncation> SmartCoreTrunc;

}

#endif // CORETRUNCATION_HPP
