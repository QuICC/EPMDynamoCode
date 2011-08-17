/** \file GenerateState.cpp
 *  \brief Initial state integrator executable
 */

// Configuration includes
//
#include "Config/Parallelisation.h"

// System includes
//
#include <iostream>

// Project includes
//
#include "General/EPMTypedefs.hpp"
#include "Domain/Truncation.hpp"
#include "Simulations/Simulation.hpp"
#include "Simulations/Traits/CodensityDiffusionTraits.hpp"
#include "Simulations/Traits/DynamoTraits.hpp"
#include "Simulations/Traits/MagneticDiffusionTraits.hpp"
#include "Simulations/Traits/ThermalRotConvectionTraits.hpp"
#include "Utilities/Integrator.hpp"

#include "Utilities/Integrals/IntegralTests.hpp"

namespace epm = EPMDynamo;

#define SIMTRAITS epm::DynamoTraits

typedef epm::IntegralTests<SIMTRAITS>  IntegralType;

typedef IntegralType::IntegralTraits  IntegralTraits;

typedef  epm::Integrator<SIMTRAITS>  Integrator;

/**
 * @brief Velocity diffusion simulation
 */
int runProgram()
{
   // Create the genertor object
   Integrator   integrator;

   // Set the codensity field on RTP decomposition
   if(SIMTRAITS::NeedCodensity && IntegralTraits::UseRTPCodensity)
   {
      IntegralType::setRTPCodensity(integrator.codC());
   }

   // Set the magetic field on RTP decomposition
   if(SIMTRAITS::NeedMagnetic && IntegralTraits::UseRTPMagnetic)
   {
      IntegralType::setRTPMagnetic(integrator.magB());
   }

   // Set the velocity field on RTP decomposition
   if(SIMTRAITS::NeedVelocity && IntegralTraits::UseRTPVelocity)
   {
      IntegralType::setRTPVelocity(integrator.velV());
   }

   // Configure transform
   integrator.configureTransform<IntegralTraits>();

   // Transform the fields
   integrator.transformRTP<IntegralTraits>();

   // Set the codensity field on Spectral decomposition
   if(SIMTRAITS::NeedCodensity && IntegralTraits::UseSpecCodensity)
   {
      IntegralType::setSpecCodensity(integrator.codC());
   }

   // Set the magetic field on Spectral decomposition
   if(SIMTRAITS::NeedMagnetic && IntegralTraits::UseSpecMagnetic)
   {
      IntegralType::setSpecMagnetic(integrator.magB());
   }

   // Set the velocity field on Spectral decomposition
   if(SIMTRAITS::NeedVelocity && IntegralTraits::UseSpecVelocity)
   {
      IntegralType::setSpecVelocity(integrator.velV());
   }

   integrator.integrate();

   // Finalise the state file
   integrator.finalise();

   return 0;
}

/**
 * @brief Main
 *
 * The actual program is in runProgram to make sure MPI initialisations
 * are called before anything else end finalization after destruction
 */
int main(int argc, char* argv[])
{
   // Initilise everything that can't be done inside a class
   epm::EPMDYNAMO_INITIALISER;

   int code;
   code = runProgram();

   // Finalise everything that can't be done inside a class
   epm::EPMDYNAMO_FINALIZER;

   return code;
}
