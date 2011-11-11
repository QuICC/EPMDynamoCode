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
#include "Simulations/Traits/VelocityDiffusionTraits.hpp"
#include "Simulations/Traits/CodensityDiffusionTraits.hpp"
#include "Simulations/Traits/MagneticDiffusionTraits.hpp"
#include "Simulations/Traits/ThermalRotConvectionTraits.hpp"
#include "Simulations/Traits/VelocityDiffusionTraits.hpp"
#include "Utilities/CSCSParaviewGenerator.hpp"

namespace epm = EPMDynamo;

#define SIMTRAITS epm::DynamoTraits

typedef  epm::CSCSParaviewGenerator<SIMTRAITS>  CSCSGenerator;

/**
 * @brief Output traits for full field values
 */
class FullVisTraits
{
   public:
      /// Output Codensity visualisation
      static const bool VisCodensity = true;

      /// Output Codensity gradient visualisation
      static const bool VisCodensityGrad = false;

      /// Output Velocity visualisation
      static const bool VisVelocity = true;

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
class TorVisTraits
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
class PolVisTraits
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
 * @brief Output traits for full field values
 */
class FullAddVisTraits
{
   public:
      /// Output Codensity visualisation
      static const bool VisCodPerturbation = false;

      /// Output Velocity visualisation
      static const bool VisInertial = true;

      /// Output Velocity only m=0 component visualisation
      static const bool VisVelM0 = true;

      /// Output Magnetic only m=0 component visualisation
      static const bool VisMagM0 = false;
};

/**
 * @brief Output traits for toroidal component values
 */
class TorAddVisTraits
{
   public:
      /// Output Codensity visualisation (doesn't make sense to be true)
      static const bool VisCodPerturbation = false;

      /// Output Velocity visualisation of toroidal component only
      static const bool VisInertial = false;

      /// Output Velocity only m=0 component visualisation
      static const bool VisVelM0 = false;

      /// Output Magnetic only m=0 component visualisation
      static const bool VisMagM0 = false;
};

/**
 * @brief Output traits for poloidal component values
 */
class PolAddVisTraits
{
   public:
      /// Output Codensity visualisation (doesn't make sense to be true)
      static const bool VisCodPerturbation = false;

      /// Output Velocity visualisation of poloidal component only
      static const bool VisInertial = false;

      /// Output Velocity only m=0 component visualisation
      static const bool VisVelM0 = false;

      /// Output Magnetic only m=0 component visualisation
      static const bool VisMagM0 = false;
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
   generator.writeTotal<FullVisTraits, FullAddVisTraits>();

   // Create visualisation output from toroidal component
   generator.writeToroidal<TorVisTraits, TorAddVisTraits>();

   // Create visualisation output from poloidal component
   generator.writePoloidal<PolVisTraits, PolAddVisTraits>();

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
