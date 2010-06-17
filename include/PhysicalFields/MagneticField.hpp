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
#include "Simulations/SimulationTraits.hpp"
#include "PhysicalFields/PhysicalVariable.hpp"
#include "PhysicalFields/PhysicalNoDivField.hpp"

namespace EPMDynamo {

   /**
    * \brief Implementation of the Magnetic vector field
    *
    * \tparam TSimType Type of the simulation
    * \tparam TFieldType Type of the field
    */
   template <typename TSimType, template <typename> class TFieldType = PhysicalNoDivField> class MagneticField: public PhysicalVariable<TSimType, TFieldType>
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

   template <typename TSimType, template <typename> class TFieldType> MagneticField<TSimType, TFieldType>::MagneticField(SmartTruncation pTrunc, typename MagneticField<TSimType, TFieldType>::TransformType &transform)
      : PhysicalVariable<TSimType, TFieldType>(pTrunc, transform)
   {
   }
}

#endif // MAGNETICFIELD_HPP
