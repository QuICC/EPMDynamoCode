/** \file SimulationConfig.hpp
 *  \brief General traits at simulation level
 */

#ifndef SIMULATIONCONFIG_HPP
#define SIMULATIONCONFIG_HPP

// System includes
//

// External includes
//

// Project includes
//

namespace EPMDynamo {

   //
   // Forward declarations

   // NumericalScheme
   class WSHSimulation;

   // Operators
   class DenseOperator;
   class DenseLUOperator;

   // TimestepScheme
   class PCSchemeTraits;

   // Parameters
   class EQRaRoParameters;
   class EEkEmRaParameters;
   class ELPmPrRaParameters;
   class EPmPrRaParameters;
   class EPmQRaParameters;
   class PrRaParameters;
   class EPmParameters;
   class EParameters;
   class InvPrRaParameters;

   // Transforms
   class MPISpectralSHTTraits;
   class SerialSpectralSHTTraits;
   class MPISHTTraits;
   class SerialSHTTraits;
   template <typename , typename > class TorPolSpectralSHTransform;

   /**
    * @brief Simulation traits class
    */
   class SimulationConfig
   {
      public:
         /// Typedef for the numerical scheme used by the solvers
         typedef WSHSimulation   NumericalScheme;

         /// Typedef for the transform type
         #ifdef EPMDYNAMO_SPLIT_RADIAL
            #ifdef EPMDYNAMO_SPLIT_SH
               typedef  TorPolSpectralSHTransform<MPISpectralSHTTraits, MPISHTTraits>   TransformType;
            #else
               typedef  TorPolSpectralSHTransform<MPISpectralSHTTraits, SerialSHTTraits>   TransformType;
            #endif //EPMDYNAMO_SPLIT_SH
         #else
            #ifdef EPMDYNAMO_SPLIT_SH
               typedef  TorPolSpectralSHTransform<SerialSpectralSHTTraits, MPISHTTraits>   TransformType;
            #else
               typedef  TorPolSpectralSHTransform<SerialSpectralSHTTraits, SerialSHTTraits>   TransformType;
            #endif //EPMDYNAMO_SPLIT_SH
         #endif //EPMDYNAMO_SPLIT_RADIAL

         /// Typedef for the operator type
         typedef  DenseOperator  OperatorType;

         /// Typedef for the factored operator type
         typedef  DenseLUOperator  FactoredOpType;

         /// Typedef for the timestepping scheme traits
         typedef  PCSchemeTraits  TimestepTraits;
//         typedef  ETDSchemeTraits  TimestepTraits;
         
         /// Typedef for the equation parameters
         //typedef EQRaRoParameters  EquationParametersType;
         typedef EPmParameters  EquationParametersType;
         //typedef EParameters  EquationParametersType;
         //typedef PrRaParameters  EquationParametersType;
         //typedef InvPrRaParameters  EquationParametersType;
   };

}

#endif // SIMULATIONCONFIG_HPP
