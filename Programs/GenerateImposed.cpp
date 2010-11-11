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
#include "Simulations/Traits/DynamoTraits.hpp"
#include "Utilities/ImposedFieldGenerator.hpp"
#include "Simulations/Types/WSHSimulation.hpp"
#include "Simulations/Types/WSHSimInc.hpp"

#include "Utilities/Imposed/BasicImposed.hpp"

namespace epm = EPMDynamo;

#define TSIMTRAITS epm::DynamoTraits
#define SIMTRAITS TSIMTRAITS<epm::WSHSimulation>

typedef epm::BasicImposed<SIMTRAITS>  ImposedType;
typedef ImposedType::ImposedTraits  ImposedTraits;

typedef  epm::ImposedFieldGenerator<epm::WSHSimulation, TSIMTRAITS>  ImposedFieldGenerator;

/**
 * @brief Velocity diffusion simulation
 */
int runProgram()
{
   // Create the genertor object
   ImposedFieldGenerator   generator;

   // Set the codensity field on RTP decomposition
   if(SIMTRAITS::NeedCodensity && ImposedTraits::UseRTPCodensity)
   {
      ImposedType::setRTPCodensity(generator.codC());
   }

   // Set the magetic field on RTP decomposition
   if(SIMTRAITS::NeedMagnetic && ImposedTraits::UseRTPMagnetic)
   {
      ImposedType::setRTPMagnetic(generator.magB());
   }

   // Set the velocity field on RTP decomposition
   if(SIMTRAITS::NeedVelocity && ImposedTraits::UseRTPVelocity)
   {
      ImposedType::setRTPVelocity(generator.velV());
   }

   // Configure transform
   generator.configureTransform<ImposedTraits>();

   // Transform the fields
   generator.transformRTP<ImposedTraits>();

   // Set the codensity field on Spectral decomposition
   if(SIMTRAITS::NeedCodensity && ImposedTraits::UseSpecCodensity)
   {
      ImposedType::setSpecCodensity(generator.codC());
   }

   // Set the magetic field on Spectral decomposition
   if(SIMTRAITS::NeedMagnetic && ImposedTraits::UseSpecMagnetic)
   {
      ImposedType::setSpecMagnetic(generator.magB());
   }

   // Set the velocity field on Spectral decomposition
   if(SIMTRAITS::NeedVelocity && ImposedTraits::UseSpecVelocity)
   {
      ImposedType::setSpecVelocity(generator.velV());
   }

   // Initialise the imposed file
   generator.setupOutput("Imposed");

   // Write the imposed file
   generator.writeImposedFile();

   // Finalise the imposed file
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
