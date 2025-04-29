/** \file CoreTruncationBase.hpp
 *  \brief Base for all the core truncation information (local or remote)
 */

#ifndef CORETRUNCATIONBASE_HPP
#define CORETRUNCATIONBASE_HPP

// System includes
//
#include <vector>

// External includes
//

// Project includes
//
#include "General/EPMTypedefs.hpp"
#include "Domain/Core/FDSHTruncation.hpp"
#include "Domain/Core/SpecTruncation.hpp"

namespace EPMDynamo {

   /**
    * \brief Implementation of the base for all core truncation information (local and remote ones)
    */
   class CoreTruncationBase
   {
      public:
         /**
          * @brief Constructor
          *
          * \param fR0 Array of radial start index
          * \param fNr Array of number of radial indexes
          * \param fM Array of harmonic orders
          * \param sL Array of harmonic degrees
          * \param sMs Vector of array of harmonic order per harmonic degree
          */
         CoreTruncationBase(const ArrayI &fR0, const ArrayI &fNr, const ArrayI &fM, const ArrayI &sL, const std::vector<ArrayI> &sMs);

         /**
          * @brief Destructor
          */
         virtual ~CoreTruncationBase() {};

         /**
          * @brief Get the RTP truncation
          */
         const SmartCoreFDSHTrunc fdsh() const;

         /**
          * @brief Get the Spec truncation
          */
         const SmartCoreSpecTrunc spec() const;
         
      protected:

      private:
         /**
          * @brief Pointer to the FDSH truncation
          */
         SmartCoreFDSHTrunc mpFDSH;

         /**
          * @brief Pointer to the Spec truncation
          */
         SmartCoreSpecTrunc mpSpec;
   };

   inline const SmartCoreFDSHTrunc CoreTruncationBase::fdsh() const
   {
      return this->mpFDSH;
   }

   inline const SmartCoreSpecTrunc CoreTruncationBase::spec() const
   {
      return this->mpSpec;
   }
}

#endif // CORETRUNCATIONBASE_HPP
