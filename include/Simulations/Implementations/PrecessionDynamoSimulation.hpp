/** \file PrecessionDynamoSimulation.hpp
 *  \brief Implementation of a dynamo simulation
 */

#ifndef PRECESSIONDYNAMOSIMULATION_HPP
#define PRECESSIONDYNAMOSIMULATION_HPP

// Configuration includes
//
#include "Config/SmartPointer.h"

// System includes
//

// External includes
//

// Project includes
//
#include "General/EPMTypedefs.hpp"
#include "Simulations/SimulationBase.hpp"
#include "Simulations/Traits/PrecessionDynamoTraits.hpp"

#include "IO/HDF5/State/StateFileReader.hpp"
#include "IO/HDF5/State/StateFileWriter.hpp"
#include "IO/HDF5/Source/SourceFileReader.hpp"
#include "IO/ASCII/EnergyFile.hpp"
#include "IO/ASCII/SpectrumFile.hpp"
#include "IO/ASCII/TimeFile.hpp"
#include "IO/ASCII/CFLTimestepFile.hpp"

#include "Equations/Induction/InductionMHD.hpp"
#include "Equations/NavierStokes/NavierStokesPrecessionMHD.hpp"

#include "BoundaryConditions/Homogeneous/ZeroBC.hpp"
#include "BoundaryConditions/Homogeneous/DDRadialBC.hpp"
#include "BoundaryConditions/Homogeneous/DRadialBC.hpp"
#include "BoundaryConditions/Homogeneous/InsulatingBC.hpp"
#include "BoundaryConditions/Homogeneous/ConductorTorBC.hpp"
#include "BoundaryConditions/Homogeneous/ConductorPolBC.hpp"
#include "BoundaryConditions/Homogeneous/StressFreeTorBC.hpp"

namespace EPMDynamo {

   /**
    * @brief Implementation of a dynamo simulation as a simulation policy class
    */
   class PrecessionDynamoSimulation: public SimulationBase
   {
      public:
         /**
          * @brief Simple empty destructor
          */
         virtual ~PrecessionDynamoSimulation() {};

      protected:
         /**
          * @brief Constructor
          */
         PrecessionDynamoSimulation();

         /**
          * @brief Initialise the fields
          */
         void initFields();

         /**
          * @brief Add HDF5 output files
          */
         void addHDF5Output();

         /**
          * @brief Add ASCII output files
          */
         void addASCIIOutput();

         /**
          * @brief Initialise the equations
          */
         void initEquations();

         /**
          * @brief Update equations RTP values
          *
          * \param step Current step in a multistep transform
          */
         void updateEquationsRTP(const int step);

         /**
          * @brief Update equations RHS
          */
         void updateEquationsRHS();

         /**
          * @brief Transform equations RHS
          *
          * \param step Current step in a multistep transform
          */
         void transformEquationsRHS(const int step);

         /**
          * @brief Include external influence
          */
         void addExternalInfluence();

         /**
          * @brief Timestep the equations
          */
         void timestepEquations();

         /**
          * @brief Configure the transforms
          */
         void configureTransforms();

      private:
         /**
          * @brief Magnetic field
          */
         PrecessionDynamoTraits::MagType   mMagB;

         /**
          * @brief Velocity field
          */
         PrecessionDynamoTraits::VelType   mVelV;

         /**
          * @brief Induction equation
          */
         InductionMHD<PrecessionDynamoTraits>    mInduction;

         /**
          * @brief Navier Stokes equation
          */
         NavierStokesPrecessionMHD<PrecessionDynamoTraits>    mNavierStokes;
   };

}

#endif // PRECESSIONDYNAMOSIMULATION_HPP
