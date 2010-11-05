/** \file GenerateCSCS.cpp
 *  \brief CSCS Paraview file generator executable
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
#include "Utilities/CSCSParaviewGenerator.hpp"
#include "Simulations/Types/WSHSimulation.hpp"
#include "Simulations/Types/WSHSimInc.hpp"

namespace epm = EPMDynamo;

#define TGENTRAITS epm::GenCBVTraits
#define GENTRAITS TGENTRAITS<epm::WSHSimulation>

typedef  epm::CSCSParaviewGenerator<epm::WSHSimulation, TGENTRAITS>  CSCSGenerator;
typedef  GENTRAITS::CodType  Codensity;
typedef  GENTRAITS::MagType  Magnetic;
typedef  GENTRAITS::VelType  Velocity;
typedef  epm::SmartTruncation  SmartTruncation;

/**
 * @brief Velocity diffusion simulation
 */
int runProgram()
{
   // Use the CSCS grid output
   epm::SimulationConstants::setCSCSGrid();

   // Create the genertor object
   CSCSGenerator   generator;

   // Initialise the state file
   generator.initOutput("CSCSParaview");

   // Transform the fields
   generator.transformSpectral();

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
