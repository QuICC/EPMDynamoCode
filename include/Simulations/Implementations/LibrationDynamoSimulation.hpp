/** \file LibrationDynamoSimulation.hpp
 *  \brief Implementation of a dynamo simulation
 */

#ifndef LIBRATIONDYNAMOSIMULATION_HPP
#define LIBRATIONDYNAMOSIMULATION_HPP

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
#include "Simulations/Traits/LibrationDynamoTraits.hpp"

#include "IO/HDF5/State/StateFileReader.hpp"
#include "IO/HDF5/State/StateFileWriter.hpp"
#include "IO/HDF5/Source/SourceFileReader.hpp"
#include "IO/ASCII/EnergyFile.hpp"
#include "IO/ASCII/SpectrumFile.hpp"
#include "IO/ASCII/TimeFile.hpp"
#include "IO/ASCII/CFLTimestepFile.hpp"

#include "Equations/Induction/InductionMHD.hpp"
#include "Equations/NavierStokes/NavierStokesLibrationMHD.hpp"

#include "BoundaryConditions/Homogeneous/ZeroBC.hpp"
#include "BoundaryConditions/Homogeneous/DDRadialBC.hpp"
#include "BoundaryConditions/Homogeneous/DRadialBC.hpp"
#include "BoundaryConditions/Homogeneous/InsulatingBC.hpp"
#include "BoundaryConditions/Homogeneous/ConductorTorBC.hpp"
#include "BoundaryConditions/Homogeneous/ConductorPolBC.hpp"
#include "BoundaryConditions/Homogeneous/StressFreeTorBC.hpp"
#include "BoundaryConditions/TimeDependent/LongitudinalLibrationBC.hpp"
#include "BoundaryConditions/TimeDependent/LatitudinalLibrationBC.hpp"

namespace EPMDynamo {

   /**
    * @brief Implementation of a dynamo simulation as a simulation policy class
    */
   class LibrationDynamoSimulation: public SimulationBase
   {
      public:
         /**
          * @brief Simple empty destructor
          */
         virtual ~LibrationDynamoSimulation() {};

      protected:
         /**
          * @brief Constructor
          */
         LibrationDynamoSimulation();

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
         LibrationDynamoTraits::MagType   mMagB;

         /**
          * @brief Velocity field
          */
         LibrationDynamoTraits::VelType   mVelV;

         /**
          * @brief Induction equation
          */
         InductionMHD<LibrationDynamoTraits>    mInduction;

         /**
          * @brief Navier Stokes equation
          */
         NavierStokesLibrationMHD<LibrationDynamoTraits>    mNavierStokes;
   };

}

#endif // LIBRATIONDYNAMOSIMULATION_HPP
