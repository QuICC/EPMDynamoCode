/** \file SimulationTraits.hpp
 *  \brief General traits at simulation level
 */

#ifndef SIMULATIONTRAITS_HPP
#define SIMULATIONTRAITS_HPP

// System includes
//
#include <boost/mpl/assert.hpp>
#include <boost/type_traits/is_same.hpp>

// External includes
//

// Project includes
//
#include "General/EPMTypedefs.hpp"
#include "Equations/Parameters/EquationParameters.hpp"
#include "Equations/Parameters/EEkEmRaParameters.hpp"
#include "Equations/Parameters/ELPmPrRaParameters.hpp"
#include "Equations/Parameters/EPmPrRaParameters.hpp"
#include "Equations/Parameters/EPmQRaParameters.hpp"
#include "Equations/Parameters/EQRaRoParameters.hpp"
#include "Timestepping/PredictorCorrector/PCSchemeTraits.hpp"
#include "Timestepping/ETD/ETDSchemeTraits.hpp"
#include "Transforms/SpectralSH/TorPolSpectralSHTransform.hpp"
#include "Transforms/SpectralSH/SerialSpectralSHTTraits.hpp"
#include "Transforms/SpectralSH/MPISpectralSHTTraits.hpp"
#include "Transforms/SphericalHarmonics/SerialSHTTraits.hpp"
#include "Transforms/SphericalHarmonics/MPISHTTraits.hpp"

namespace EPMDynamo {

   // Forward declaration of the Operators
   class DenseOperator;
   class DenseLUOperator;

   /**
    * @brief Simulation traits class
    *
    * \tparam TSimType Type of the simulation
    */
   template <typename TSimType> class SimulationTraits
   {
      public:
         /// Typedef for the transform type
         #ifdef EPMDYNAMO_SPLIT_RADIAL
            #ifdef EPMDYNAMO_SPLIT_SH
               typedef  TorPolSpectralSHTransform<TSimType, MPISpectralSHTTraits, MPISHTTraits>   TransformType;
            #else
               typedef  TorPolSpectralSHTransform<TSimType, MPISpectralSHTTraits, SerialSHTTraits>   TransformType;
            #endif //EPMDYNAMO_SPLIT_SH
         #else
            #ifdef EPMDYNAMO_SPLIT_SH
               typedef  TorPolSpectralSHTransform<TSimType, SerialSpectralSHTTraits, MPISHTTraits>   TransformType;
            #else
               typedef  TorPolSpectralSHTransform<TSimType, SerialSpectralSHTTraits, SerialSHTTraits>   TransformType;
            #endif //EPMDYNAMO_SPLIT_SH
         #endif //EPMDYNAMO_SPLIT_RADIAL

         /// Typedef for the operator type
         typedef  DenseOperator  OperatorType;

         /// Typedef for the factored operator type
         typedef  DenseLUOperator  FactoredOpType;

         /// Typedef for the timestepping scheme traits
//         typedef  PCSchemeTraits<TSimType>  TimestepTraits;
         typedef  ETDSchemeTraits<TSimType>  TimestepTraits;
         
         /// Typedef for the equation parameters
         typedef EQRaRoParameters  EquationParametersType;

         // Make sure the non finished parameters implementation are not used: EEkEmRa case
         BOOST_MPL_ASSERT_NOT(( boost::is_same< EquationParametersType,EEkEmRaParameters > ));
         // Make sure the non finished parameters implementation are not used: ELPmPrRa case
         BOOST_MPL_ASSERT_NOT(( boost::is_same< EquationParametersType,ELPmPrRaParameters > ));
         // Make sure the non finished parameters implementation are not used: EPmPrRa case
         BOOST_MPL_ASSERT_NOT(( boost::is_same< EquationParametersType,EPmPrRaParameters > ));
         // Make sure the non finished parameters implementation are not used: EPmQRa case
         BOOST_MPL_ASSERT_NOT(( boost::is_same< EquationParametersType,EPmQRaParameters > ));
   };

}

#endif // SIMULATIONTRAITS_HPP
