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

#define TGENTRAITS epm::GenCBVTraits
#define GENTRAITS TGENTRAITS<epm::WSHSimulation>

//typedef epm::RandomState<GENTRAITS>  StateType;
typedef epm::SakurabaState<GENTRAITS>  StateType;

typedef  epm::InitialStateGenerator<epm::WSHSimulation, TGENTRAITS>  IStateGenerator;
typedef  GENTRAITS::CodType  Codensity;
typedef  GENTRAITS::MagType  Magnetic;
typedef  GENTRAITS::VelType  Velocity;
typedef  epm::SmartTruncation  SmartTruncation;

/**
 * @brief Velocity diffusion simulation
 */
int runProgram()
{
   // Create the genertor object
   IStateGenerator   generator;

   // Set the codensity field
   if(GENTRAITS::NeedCodensity)
   {
      StateType::setRTPCodensity(generator.codC());
   }

   // Set the magetic field
   if(GENTRAITS::NeedMagnetic)
   {
      StateType::setRTPMagnetic(generator.magB());
   }

   // Set the velocity field
   if(GENTRAITS::NeedVelocity)
   {
      StateType::setRTPVelocity(generator.velV());
   }

   // Transform the fields
   generator.transformRTP();

   // Set the codensity field
   if(GENTRAITS::NeedCodensity)
   {
      StateType::setSpecCodensity(generator.codC());
   }

   // Set the magetic field
   if(GENTRAITS::NeedMagnetic)
   {
      StateType::setSpecMagnetic(generator.magB());
   }

   // Set the velocity field
   if(GENTRAITS::NeedVelocity)
   {
      StateType::setSpecVelocity(generator.velV());
   }

   // Initialise the state file
   generator.initOutput("Initial");

   // Write the state file
   generator.writeFile();

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
