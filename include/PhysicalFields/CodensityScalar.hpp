/** \file CodensityScalar.hpp
 *  \brief Implementation of the Codensity scalar field
 */

#ifndef CODENSITYSCALAR_HPP
#define CODENSITYSCALAR_HPP

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
#include "PhysicalFields/Traits/DefaultCodTraits.hpp"
#include "PhysicalFields/Traits/ScalarNTermsTraits.hpp"

namespace EPMDynamo {

   /**
    * \brief Implementation of the Codensity scalar
    *
    * \tparam TCodTraits Type of the field
    */
   template <typename TCodTraits = DefaultCodTraits> class CodensityScalar: public PhysicalVariable<TCodTraits>
   {
      public:
         /// Typedef from Simulation trait to local transform type
         typedef SimulationConfig::TransformType    TransformType;

         /// Typedef for the type of the required non-linear terms
         typedef PhysicalVariable<ScalarNTermsTraits>    NTermsType;

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


   template <typename TCodTraits> CodensityScalar<TCodTraits>::CodensityScalar(SmartTruncation pTrunc, typename CodensityScalar<TCodTraits>::TransformType &transform)
      : PhysicalVariable<TCodTraits>(pTrunc, transform)
   {
   }

}

#endif // CODENSITYSCALAR_HPP
