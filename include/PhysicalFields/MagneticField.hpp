/** \file MagneticField.hpp
 *  \brief Implementation of the Magnetic vector field
 */

#ifndef MAGNETICFIELD_HPP
#define MAGNETICFIELD_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "Simulations/Traits/SimulationTraits.hpp"
#include "PhysicalFields/PhysicalVariable.hpp"
#include "PhysicalFields/Traits/DefaultMagTraits.hpp"

namespace EPMDynamo {

   /**
    * \brief Implementation of the Magnetic vector field
    *
    * \tparam TSimType Type of the simulation
    * \tparam TMagTraits Type of the field
    */
   template <typename TSimType, template <typename> class TMagTraits = DefaultMagTraits> class MagneticField: public PhysicalVariable<TSimType, TMagTraits>
   {
      public:
         /// Typedef from Simulation trait to local transform type
         typedef typename SimulationTraits<TSimType>::TransformType    TransformType;

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

   template <typename TSimType, template <typename> class TMagTraits> MagneticField<TSimType, TMagTraits>::MagneticField(SmartTruncation pTrunc, typename MagneticField<TSimType, TMagTraits>::TransformType &transform)
      : PhysicalVariable<TSimType, TMagTraits>(pTrunc, transform)
   {
   }
}

#endif // MAGNETICFIELD_HPP
