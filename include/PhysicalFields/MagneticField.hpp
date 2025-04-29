/** \file MagneticField.hpp
 *  \brief Implementation of the Magnetic vector field
 */

#ifndef MAGNETICFIELD_HPP
#define MAGNETICFIELD_HPP

// Configuration includes
//
#include "Config/SimulationConfig.hpp"

// System includes
//

// External includes
//

// Project includes
//
#include "PhysicalFields/PhysicalVariable.hpp"
#include "PhysicalFields/Traits/DefaultMagTraits.hpp"
#include "PhysicalFields/Traits/FieldNTermsTraits.hpp"

namespace EPMDynamo {

   /**
    * \brief Implementation of the Magnetic vector field
    *
    * \tparam TMagTraits Type of the field
    */
   template <typename TMagTraits = DefaultMagTraits> class MagneticField: public PhysicalVariable<TMagTraits>
   {
      public:
         /// Typedef from Simulation trait to local transform type
         typedef SimulationConfig::TransformType    TransformType;

         /// Typedef for the type of the required non-linear terms
         typedef PhysicalVariable<FieldNTermsTraits>    NTermsType;

         /**
          * @brief Constructor
          *
          * @param pTrunc Truncation information
          * @param transform Reference to the transform object
          */
         MagneticField(SmartTruncation pTrunc, TransformType &transform);

         /**
          * @brief Simple empty destructor
          */
         virtual ~MagneticField() {};

      protected:

      private:
   };

   template <typename TMagTraits> MagneticField<TMagTraits>::MagneticField(SmartTruncation pTrunc, typename MagneticField<TMagTraits>::TransformType &transform)
      : PhysicalVariable<TMagTraits>(pTrunc, transform)
   {
   }
}

#endif // MAGNETICFIELD_HPP
