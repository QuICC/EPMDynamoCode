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
#include "Simulations/Traits/DynamoTraits.hpp"
#include "Simulations/Traits/ThermalConvectionTraits.hpp"
#include "Simulations/Traits/ThermalRotConvectionTraits.hpp"
#include "Simulations/Traits/ThermalNLRotConvectionTraits.hpp"

#include "Utilities/SourceGenerator.hpp"

#include "Utilities/Sources/JonesSource.hpp"
#include "Utilities/Sources/SakurabaSource.hpp"
#include "Utilities/Sources/HomogeneousSource.hpp"

namespace epm = EPMDynamo;


//#define SIMTRAITS epm::ThermalConvectionTraits
//#define SIMTRAITS epm::ThermalRotConvectionTraits
#define SIMTRAITS epm::ThermalNLRotConvectionTraits
//#define SIMTRAITS epm::VelocityNLRotDiffusionTraits

typedef epm::HomogeneousSource<SIMTRAITS>  SourceType;
typedef SourceType::SourceTraits  SourceTraits;

typedef  epm::SourceGenerator<SIMTRAITS>  SourceGenerator;

/**
 * @brief Velocity diffusion simulation
 */
int runProgram()
{
   // Create the genertor object
   SourceGenerator   generator;

   // Set the codensity field on RTP decomposition
   if(SIMTRAITS::NeedCodensity && SourceTraits::UseRTPCodensity)
   {
      SourceType::setRTPCodensity(generator.codC());
   }

   // Set the magetic field on RTP decomposition
   if(SIMTRAITS::NeedMagnetic && SourceTraits::UseRTPMagnetic)
   {
      SourceType::setRTPMagnetic(generator.magB());
   }

   // Set the velocity field on RTP decomposition
   if(SIMTRAITS::NeedVelocity && SourceTraits::UseRTPVelocity)
   {
      SourceType::setRTPVelocity(generator.velV());
   }

   // Configure transform
   generator.configureTransform<SourceTraits>();

   // Transform the fields
   generator.transformRTP<SourceTraits>();

   // Set the codensity field on Spectral decomposition
   if(SIMTRAITS::NeedCodensity && SourceTraits::UseSpecCodensity)
   {
      SourceType::setSpecCodensity(generator.codC());
   }

   // Set the magetic field on Spectral decomposition
   if(SIMTRAITS::NeedMagnetic && SourceTraits::UseSpecMagnetic)
   {
      SourceType::setSpecMagnetic(generator.magB());
   }

   // Set the velocity field on Spectral decomposition
   if(SIMTRAITS::NeedVelocity && SourceTraits::UseSpecVelocity)
   {
      SourceType::setSpecVelocity(generator.velV());
   }

   // Initialise the source file
   generator.setupOutput("Source");

   // Write the source file
   generator.writeSourceFile();

   // Finalise the source file
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
