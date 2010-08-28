/** \file GenerateSource.cpp
 *  \brief Source file generator executable
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
#include "Utilities/SourceGenerator.hpp"
#include "Simulations/Types/WSHSimulation.hpp"
#include "Simulations/Types/WSHSimInc.hpp"

#include "Utilities/Sources/SakurabaSource.hpp"

namespace epm = EPMDynamo;

#define TGENTRAITS epm::GenCTraits
#define GENTRAITS TGENTRAITS<epm::WSHSimulation>

typedef epm::SakurabaSource<GENTRAITS>  SourceType;

typedef  epm::SourceGenerator<epm::WSHSimulation, TGENTRAITS>  SourceGenerator;
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
   SourceGenerator   generator;

   // Set the codensity field on RTP decomposition
   if(GENTRAITS::NeedCodensity)
   {
      SourceType::setRTPCodensity(generator.codC());
   }

   // Set the magetic field on RTP decomposition
   if(GENTRAITS::NeedMagnetic)
   {
      SourceType::setRTPMagnetic(generator.magB());
   }

   // Set the velocity field on RTP decomposition
   if(GENTRAITS::NeedVelocity)
   {
      SourceType::setRTPVelocity(generator.velV());
   }

   // Transform the fields
   generator.transformRTP();

   // Set the codensity field on Spectral decomposition
   if(GENTRAITS::NeedCodensity)
   {
      SourceType::setSpecCodensity(generator.codC());
   }

   // Set the magetic field on Spectral decomposition
   if(GENTRAITS::NeedMagnetic)
   {
      SourceType::setSpecMagnetic(generator.magB());
   }

   // Set the velocity field on Spectral decomposition
   if(GENTRAITS::NeedVelocity)
   {
      SourceType::setSpecVelocity(generator.velV());
   }

   // Initialise the state file
   generator.initOutput("Source");

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
