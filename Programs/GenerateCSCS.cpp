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
#include "Simulations/Traits/DynamoTraits.hpp"
#include "Simulations/Traits/MagneticDiffusionTraits.hpp"
#include "Simulations/Traits/ThermalRotConvectionTraits.hpp"
#include "Utilities/CSCSParaviewGenerator.hpp"

namespace epm = EPMDynamo;

#define SIMTRAITS epm::MagneticDiffusionTraits

typedef  epm::CSCSParaviewGenerator<SIMTRAITS>  CSCSGenerator;

/**
 * @brief Output traits for full field values
 */
class FullFieldTraits
{
   public:
      /// Output Codensity visualisation
      static const bool VisCodensity = false;

      /// Output Codensity gradient visualisation
      static const bool VisCodensityGrad = false;

      /// Output Velocity visualisation
      static const bool VisVelocity = false;

      /// Output Vorticity visualisation
      static const bool VisVorticity = false;

      /// Output Helicity visualisation
      static const bool VisHelicity = false;

      /// Output Magnetic visualisation
      static const bool VisMagnetic = true;

      /// Output Lorentz visualisation
      static const bool VisLorentz = true;
};

/**
 * @brief Output traits for toroidal component values
 */
class ToroidalTraits
{
   public:
      /// Output Codensity visualisation (doesn't make sense to be true)
      static const bool VisCodensity = false;

      /// Output Codensity gradient visualisation (doesn't make sense to be true)
      static const bool VisCodensityGrad = false;

      /// Output Velocity visualisation of toroidal component only
      static const bool VisVelocity = false;

      /// Output Vorticity visualisation of toroidal component only
      static const bool VisVorticity = false;

      /// Output Helicity visualisation of toroidal component only
      static const bool VisHelicity = false;

      /// Output Magnetic visualisation of toroidal component only
      static const bool VisMagnetic = true;

      /// Output Lorentz visualisation of toroidal component only
      static const bool VisLorentz = false;
};

/**
 * @brief Output traits for poloidal component values
 */
class PoloidalTraits
{
   public:
      /// Output Codensity visualisation (doesn't make sense to be true)
      static const bool VisCodensity = false;

      /// Output Codensity gradient visualisation (doesn't make sense to be true)
      static const bool VisCodensityGrad = false;

      /// Output Velocity visualisation of poloidal component only
      static const bool VisVelocity = false;

      /// Output Vorticity visualisation of poloidal component only
      static const bool VisVorticity = false;

      /// Output Helicity visualisation of poloidal component only
      static const bool VisHelicity = false;

      /// Output Magnetic visualisation of poloidal component only
      static const bool VisMagnetic = true;

      /// Output Lorentz visualisation of poloidal component only
      static const bool VisLorentz = false;
};

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
   generator.setupOutput("CSCSParaview");

   // Create visualisation output from total field
   generator.writeTotal<FullFieldTraits>();

   // Create visualisation output from toroidal component
   generator.writeToroidal<ToroidalTraits>();

   // Create visualisation output from poloidal component
   generator.writePoloidal<PoloidalTraits>();

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
