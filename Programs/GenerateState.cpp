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

namespace epm = EPMDynamo;

#define TGENTRAITS epm::GenCBVTraits
#define GENTRAITS TGENTRAITS<epm::WSHSimulation>

typedef  epm::InitialStateGenerator<epm::WSHSimulation, TGENTRAITS>  IStateGenerator;
typedef  GENTRAITS::CodType  Codensity;
typedef  GENTRAITS::MagType  Magnetic;
typedef  GENTRAITS::VelType  Velocity;
typedef  epm::SmartTruncation  SmartTruncation;

void setRTPCodensity(Codensity &codC)
{
   SmartTruncation pTrunc = codC.oc().trunc();

   for(int n=0; n < pTrunc->local()->rtp()->nR(); ++n)
   {
      codC.rOc().rRTP().rShell(n).setConstant(0.0);
   }
}

void setRTPMagnetic(Magnetic &magB)
{
   SmartTruncation pTrunc = magB.oc().trunc();

   for(int n=0; n < pTrunc->local()->rtp()->nR(); ++n)
   {
      magB.rOc().rRTP().rR().rShell(n).setConstant(0.0);
      magB.rOc().rRTP().rTheta().rShell(n).setConstant(0.0);
      magB.rOc().rRTP().rPhi().rShell(n).setConstant(0.0);
   }
}

void setRTPVelocity(Velocity &velV)
{
   SmartTruncation pTrunc = velV.oc().trunc();

   for(int n=0; n < pTrunc->local()->rtp()->nR(); ++n)
   {
      velV.rOc().rRTP().rR().rShell(n).setConstant(0.);
      velV.rOc().rRTP().rTheta().rShell(n).setConstant(0.);
      velV.rOc().rRTP().rPhi().rShell(n).setConstant(0.);
   }
}

void setSpecCodensity(Codensity &codC)
{
   SmartTruncation pTrunc = codC.oc().trunc();

   codC.rOc().rPerturbation().rLShell(0)(0,0) = 1.0/4.0;
   codC.rOc().rPerturbation().rLShell(0)(1,0) = -1.0/2.0;
}

void setSpecMagnetic(Magnetic &magB)
{
   SmartTruncation pTrunc = magB.oc().trunc();

   for(int l=1; l < pTrunc->local()->spec()->nL()/2; ++l)
   {
      // Set some perturbation random energy in Toroidal component
      magB.rOc().rPerturbation().rTor().rLShell(l).block(0,0,pTrunc->sim()->rad()->nN()/2, pTrunc->local()->spec()->nM(l)/2).setRandom();
      magB.rOc().rPerturbation().rTor().rLShell(l) *= epm::EPMComplex(1.0e-5,0.0);

      // Set some perturbation random energy in Poloidal component
      magB.rOc().rPerturbation().rPol().rLShell(l).block(0,0,pTrunc->sim()->rad()->nN()/2, pTrunc->local()->spec()->nM(l)/2).setRandom();
      magB.rOc().rPerturbation().rPol().rLShell(l) *= epm::EPMComplex(1.0e-5,0.0);
   }
}

void setSpecVelocity(Velocity &velV)
{
   SmartTruncation pTrunc = velV.oc().trunc();

   for(int l=1; l < pTrunc->local()->spec()->nL()/2; ++l)
   {
      // Set some perturbation random energy in Toroidal component
      velV.rOc().rPerturbation().rTor().rLShell(l).block(0,0,pTrunc->sim()->rad()->nN()/2, pTrunc->local()->spec()->nM(l)/2).setRandom();
      velV.rOc().rPerturbation().rTor().rLShell(l) *= epm::EPMComplex(1.0e-5,0.0);

      // Set some perturbation random energy in Poloidal component
      velV.rOc().rPerturbation().rPol().rLShell(l).block(0,0,pTrunc->sim()->rad()->nN()/2, pTrunc->local()->spec()->nM(l)/2).setRandom();
      velV.rOc().rPerturbation().rPol().rLShell(l) *= epm::EPMComplex(1.0e-5, 0.0);
   }
}

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

   // Set the codensity field
   if(GENTRAITS::NeedCodensity)
   {
      setSpecCodensity(generator.codC());
   }

   // Set the magetic field
   if(GENTRAITS::NeedMagnetic)
   {
      setSpecMagnetic(generator.magB());
   }

   // Set the velocity field
   if(GENTRAITS::NeedVelocity)
   {
      setSpecVelocity(generator.velV());
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
