/** \file AnalyseState.cpp
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
#include "Simulations/Traits/DynamoTraits.hpp"
#include "Simulations/Traits/CodensityDiffusionTraits.hpp"
#include "Simulations/Traits/MagneticDiffusionTraits.hpp"
#include "Simulations/Traits/ThermalRotConvectionTraits.hpp"
#include "Simulations/Traits/VelocityDiffusionTraits.hpp"
#include "Utilities/StateAnalyser.hpp"

#include "Utilities/Analysis/LibrationAnalysis.hpp"
#include "Utilities/Analysis/ExtractRadiusAnalysis.hpp"

namespace epm = EPMDynamo;

#define SIMTRAITS epm::VelocityDiffusionTraits

typedef  epm::StateAnalyser<SIMTRAITS, epm::ExtractRadiusAnalysis>  Analyser;

/**
 * @brief Velocity diffusion simulation
 */
int runProgram()
{
   // Create the analyser object
   Analyser   analyser;

   // Initialise the analyser
   analyser.setup("2Analyse");

   // Configure transform
   analyser.configureTransform<Analyser::AnalysisTraits>();

   // Transform the fields
   analyser.transformSpectral<Analyser::AnalysisTraits>();

   // Analyse the codensity field RTP decomposition
   if(SIMTRAITS::NeedCodensity && Analyser::AnalysisTraits::UseRTPCodensity)
   {
      analyser.analyseRTPCodensity();
   }

   // Set the magetic field on RTP decomposition
   if(SIMTRAITS::NeedMagnetic && Analyser::AnalysisTraits::UseRTPMagnetic)
   {
      analyser.analyseRTPMagnetic();
   }

   // Set the velocity field on RTP decomposition
   if(SIMTRAITS::NeedVelocity && Analyser::AnalysisTraits::UseRTPVelocity)
   {
      analyser.analyseRTPVelocity();
   }

   // Set the codensity field on Spectral decomposition
   if(SIMTRAITS::NeedCodensity && Analyser::AnalysisTraits::UseSpecCodensity)
   {
      analyser.analyseSpecCodensity();
   }

   // Set the magetic field on Spectral decomposition
   if(SIMTRAITS::NeedMagnetic && Analyser::AnalysisTraits::UseSpecMagnetic)
   {
      analyser.analyseSpecMagnetic();
   }

   // Set the velocity field on Spectral decomposition
   if(SIMTRAITS::NeedVelocity && Analyser::AnalysisTraits::UseSpecVelocity)
   {
      analyser.analyseSpecVelocity();
   }

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
