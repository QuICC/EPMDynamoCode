/** \file CodensityScalar.hpp
 *  \brief Implementation of the Codensity scalar field
 */

#ifndef CODENSITYSCALAR_HPP
#define CODENSITYSCALAR_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "Simulations/Traits/SimulationTraits.hpp"
#include "PhysicalFields/PhysicalScalar.hpp"
#include "PhysicalFields/PhysicalVariable.hpp"
#include "PhysicalFields/Extensions/WithSpectralSource.hpp"

namespace EPMDynamo {

   /**
    * \brief Implementation of the Codensity scalar
    *
    * \tparam TSimType Type of the simulation
    * \tparam TFieldType Type of the field
    */
   template <typename TSimType, template <typename> class TFieldType = PhysicalScalar> class CodensityScalar: public PhysicalVariable<TSimType, TFieldType, WithSpectralSource>
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
         CodensityScalar(SmartTruncation pTrunc, TransformType &transform);

         /**
          * @brief Destructor
          */
         virtual ~CodensityScalar() {};
         
      protected:

      private:

   };


   template <typename TSimType, template <typename> class TFieldType> CodensityScalar<TSimType, TFieldType>::CodensityScalar(SmartTruncation pTrunc, typename CodensityScalar<TSimType, TFieldType>::TransformType &transform)
      : PhysicalVariable<TSimType, TFieldType, WithSpectralSource>(pTrunc, transform)
   {
   }

}

#endif // CODENSITYSCALAR_HPP
