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

namespace epm = EPMDynamo;

#define TGENTRAITS epm::GenCVTraits
#define GENTRAITS TGENTRAITS<epm::WSHSimulation>

typedef  epm::ImposedFieldGenerator<epm::WSHSimulation, TGENTRAITS>  ImposedFieldGenerator;
typedef  GENTRAITS::CodType  Codensity;
typedef  GENTRAITS::MagType  Magnetic;
typedef  GENTRAITS::VelType  Velocity;
typedef  epm::SmartTruncation  SmartTruncation;

void setRTPCodensity(Codensity &codC)
{
   SmartTruncation pTrunc = codC.trunc();

   for(int n=0; n < pTrunc->local()->rtp()->nR(); ++n)
   {
      codC.rOc().rRTP().rShell(n).setConstant(0.0);
   }
}

void setRTPMagnetic(Magnetic &magB)
{
   SmartTruncation pTrunc = magB.trunc();

   for(int n=0; n < pTrunc->local()->rtp()->nR(); ++n)
   {
      magB.rOc().rRTP().rR().rShell(n).setConstant(0.0);
      magB.rOc().rRTP().rTheta().rShell(n).setConstant(0.0);
      magB.rOc().rRTP().rPhi().rShell(n).setConstant(0.0);
   }
}

void setRTPVelocity(Velocity &velV)
{
   SmartTruncation pTrunc = velV.trunc();

   for(int n=0; n < pTrunc->local()->rtp()->nR(); ++n)
   {
      velV.rOc().rRTP().rR().rShell(n).setConstant(0.);
      velV.rOc().rRTP().rTheta().rShell(n).setConstant(0.);
      velV.rOc().rRTP().rPhi().rShell(n).setConstant(0.);
   }
}

/**
 * @brief Velocity diffusion simulation
 */
int runProgram()
{
   // Create the genertor object
   ImposedFieldGenerator   generator;

   // Set the codensity field
   if(GENTRAITS::NeedCodensity)
   {
      setRTPCodensity(generator.codC());
   }

   // Set the magetic field
   if(GENTRAITS::NeedMagnetic)
   {
      setRTPMagnetic(generator.magB());
   }

   // Set the velocity field
   if(GENTRAITS::NeedVelocity)
   {
      setRTPVelocity(generator.velV());
   }

   // Transform the fields
   generator.transformRTP();

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
