/** \file FieldNTermsFieldTraits.hpp
 *  \brief Traits for the field of the vector field non linear terms
 */

#ifndef FIELDNTERMSFIELDTRAITS_HPP
#define FIELDNTERMSFIELDTRAITS_HPP

// Configuration includes
//

// System includes
//

// External includes
//

// Project includes
//
#include "PhysicalFields/Fields/PhysicalTorPolBase.hpp"
#include "PhysicalFields/Fields/PhysicalTorPolField.hpp"

namespace EPMDynamo {

   /**
    * \brief Traits for the field of the vector field non linear terms
    */
   class FieldNTermsFieldTraits
   {
      public:
         /// Typedef for the type of the spectral field
         typedef PhysicalTorPolField<PhysicalTorPolBase>::SpectralFieldType  SpectralFieldType;

         /// Typedef for the type of an inner core field
         typedef PhysicalTorPolField<PhysicalTorPolBase> ICFieldType;

         /// Typedef for the type of an outer core field
         typedef PhysicalTorPolField<PhysicalTorPolBase> OCFieldType;
   };

}

#endif // FIELDNTERMSFIELDTRAITS_HPP
