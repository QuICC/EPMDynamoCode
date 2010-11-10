/** \file GenerateState.cpp
 *  \brief Initial state generator executable
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
#include "Simulations/Traits/DynamoTraits.hpp"
#include "Utilities/Traits/GenCBVTraits.hpp"
#include "Utilities/Traits/GenCVTraits.hpp"
#include "Utilities/Traits/GenBVTraits.hpp"
#include "Utilities/Traits/GenCTraits.hpp"
#include "Utilities/Traits/GenBTraits.hpp"
#include "Utilities/Traits/GenVTraits.hpp"
#include "Utilities/InitialStateGenerator.hpp"
#include "Simulations/Types/WSHSimulation.hpp"
#include "Simulations/Types/WSHSimInc.hpp"

#include "Utilities/States/RandomState.hpp"
#include "Utilities/States/SakurabaState.hpp"

namespace epm = EPMDynamo;

#define TSIMTRAITS epm::DynamoTraits
#define SIMTRAITS TSIMTRAITS<epm::WSHSimulation>

//typedef epm::RandomState<SIMTRAITS>  StateType;
typedef epm::SakurabaState<SIMTRAITS>  StateType;

typedef StateType::StateTraits  StateTraits;

typedef  epm::InitialStateGenerator<epm::WSHSimulation, TSIMTRAITS>  IStateGenerator;

/**
 * @brief Velocity diffusion simulation
 */
int runProgram()
{
   // Create the genertor object
   IStateGenerator   generator;

   // Set the codensity field on RTP decomposition
   if(SIMTRAITS::NeedCodensity && StateTraits::UseRTPCodensity)
   {
      StateType::setRTPCodensity(generator.codC());
   }

   // Set the magetic field on RTP decomposition
   if(SIMTRAITS::NeedMagnetic && StateTraits::UseRTPMagnetic)
   {
      StateType::setRTPMagnetic(generator.magB());
   }

   // Set the velocity field on RTP decomposition
   if(SIMTRAITS::NeedVelocity && StateTraits::UseRTPVelocity)
   {
      StateType::setRTPVelocity(generator.velV());
   }

   // Configure transform
   generator.configureTransform<StateTraits>();

   // Transform the fields
   generator.transformRTP<StateTraits>();

   // Set the codensity field on Spectral decomposition
   if(SIMTRAITS::NeedCodensity && StateTraits::UseSpecCodensity)
   {
      StateType::setSpecCodensity(generator.codC());
   }

   // Set the magetic field on Spectral decomposition
   if(SIMTRAITS::NeedMagnetic && StateTraits::UseSpecMagnetic)
   {
      StateType::setSpecMagnetic(generator.magB());
   }

   // Set the velocity field on Spectral decomposition
   if(SIMTRAITS::NeedVelocity && StateTraits::UseSpecVelocity)
   {
      StateType::setSpecVelocity(generator.velV());
   }

   // Initialise the state file
   generator.setupOutput("Initial");

   // Write the state file
   generator.writeStateFile();

   // Finalise the state file
   generator.finalise();

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
