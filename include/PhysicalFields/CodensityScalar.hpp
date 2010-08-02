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
#include "PhysicalFields/PhysicalVariable.hpp"
#include "PhysicalFields/Traits/DefaultCodTraits.hpp"

namespace EPMDynamo {

   /**
    * \brief Implementation of the Codensity scalar
    *
    * \tparam TSimType Type of the simulation
    * \tparam TCodTraits Type of the field
    */
   template <typename TSimType, template <typename> class TCodTraits = DefaultCodTraits> class CodensityScalar: public PhysicalVariable<TSimType, TCodTraits>
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


   template <typename TSimType, template <typename> class TCodTraits> CodensityScalar<TSimType, TCodTraits>::CodensityScalar(SmartTruncation pTrunc, typename CodensityScalar<TSimType, TCodTraits>::TransformType &transform)
      : PhysicalVariable<TSimType, TCodTraits>(pTrunc, transform)
   {
   }

}

#endif // CODENSITYSCALAR_HPP
