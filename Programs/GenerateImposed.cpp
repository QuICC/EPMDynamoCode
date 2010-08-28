/** \file GenerateImposed.cpp
 *  \brief Imposed field file generator executable
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
#include "Utilities/ImposedFieldGenerator.hpp"
#include "Simulations/Types/WSHSimulation.hpp"
#include "Simulations/Types/WSHSimInc.hpp"

#include "Utilities/Imposed/BasicImposed.hpp"

namespace epm = EPMDynamo;

#define TGENTRAITS epm::GenCVTraits
#define GENTRAITS TGENTRAITS<epm::WSHSimulation>

typedef epm::BasicImposed<GENTRAITS>  ImposedType;

typedef  epm::ImposedFieldGenerator<epm::WSHSimulation, TGENTRAITS>  ImposedFieldGenerator;
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
   ImposedFieldGenerator   generator;

   // Set the codensity field on RTP decomposition
   if(GENTRAITS::NeedCodensity)
   {
      ImposedType::setRTPCodensity(generator.codC());
   }

   // Set the magetic field on RTP decomposition
   if(GENTRAITS::NeedMagnetic)
   {
      ImposedType::setRTPMagnetic(generator.magB());
   }

   // Set the velocity field on RTP decomposition
   if(GENTRAITS::NeedVelocity)
   {
      ImposedType::setRTPVelocity(generator.velV());
   }

   // Transform the fields
   generator.transformRTP();

   // Set the codensity field on Spectral decomposition
   if(GENTRAITS::NeedCodensity)
   {
      ImposedType::setSpecCodensity(generator.codC());
   }

   // Set the magetic field on Spectral decomposition
   if(GENTRAITS::NeedMagnetic)
   {
      ImposedType::setSpecMagnetic(generator.magB());
   }

   // Set the velocity field on Spectral decomposition
   if(GENTRAITS::NeedVelocity)
   {
      ImposedType::setSpecVelocity(generator.velV());
   }

   // Initialise the state file
   generator.initOutput("Imposed");

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
