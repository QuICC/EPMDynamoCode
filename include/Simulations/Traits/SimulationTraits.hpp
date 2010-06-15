/** \file SimulationTraits.hpp
 *  \brief General traits at simulation level
 */

#ifndef SIMULATIONTRAITS_HPP
#define SIMULATIONTRAITS_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "General/Typedefs.hpp"
#include "Timestepping/PredictorCorrector/PCSchemeTraits.hpp"
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
         #ifdef EPMDYNAMO_RADIAL_MPI
            #ifdef EPMDYNAMO_SH_MPI
               typedef  TorPolSpectralSHTransform<TSimType, MPISpectralSHTTraits, MPISHTTraits>   TransformType;
            #else
               typedef  TorPolSpectralSHTransform<TSimType, MPISpectralSHTTraits, SerialSHTTraits>   TransformType;
            #endif //EPMDYNAMO_SH_MPI
         #else
            #ifdef EPMDYNAMO_SH_MPI
               typedef  TorPolSpectralSHTransform<TSimType, SerialSpectralSHTTraits, MPISHTTraits>   TransformType;
            #else
               typedef  TorPolSpectralSHTransform<TSimType, SerialSpectralSHTTraits, SerialSHTTraits>   TransformType;
            #endif //EPMDYNAMO_SH_MPI
         #endif //EPMDYNAMO_RADIAL_MPI

         /// Typedef for the operator type
         typedef  DenseOperator  OperatorType;

         /// Typedef for the factored operator type
         typedef  DenseLUOperator  FactoredOpType;

         /// Typedef for the timestepping scheme traits
         typedef  PCSchemeTraits<TSimType>  TimestepTraits;
   };

}

#endif // SIMULATIONTRAITS_HPP
