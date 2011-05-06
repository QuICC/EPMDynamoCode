/** \file DynamoSimulation.hpp
 *  \brief Implementation of a dynamo simulation
 */

#ifndef DYNAMOSIMULATION_HPP
#define DYNAMOSIMULATION_HPP

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
#include "Simulations/Traits/DynamoTraits.hpp"

#include "IO/HDF5/State/StateFileReader.hpp"
#include "IO/HDF5/State/StateFileWriter.hpp"
#include "IO/HDF5/Source/SourceFileReader.hpp"
#include "IO/ASCII/EnergyFile.hpp"
#include "IO/ASCII/SpectrumFile.hpp"
#include "IO/ASCII/TimeFile.hpp"
#include "IO/ASCII/CFLTimestepFile.hpp"

#include "Equations/Induction/InductionMHD.hpp"
#include "Equations/Transport/TransportMHD.hpp"
#include "Equations/NavierStokes/NavierStokesMHD.hpp"

#include "BoundaryConditions/Homogeneous/ZeroBC.hpp"
#include "BoundaryConditions/Homogeneous/DDRadialBC.hpp"
#include "BoundaryConditions/Homogeneous/DRadialBC.hpp"
#include "BoundaryConditions/Homogeneous/InsulatingBC.hpp"
#include "BoundaryConditions/L0Harmonic/DifferentialHeatingBC.hpp"
#include "BoundaryConditions/Homogeneous/StressFreeTorBC.hpp"

namespace EPMDynamo {

   /**
    * @brief Implementation of a dynamo simulation as a simulation policy class
    */
   class DynamoSimulation: public SimulationBase
   {
      public:
         /**
          * @brief Simple empty destructor
          */
         virtual ~DynamoSimulation() {};

      protected:
         /**
          * @brief Constructor
          */
         DynamoSimulation();

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
          * @brief Codensity scalar
          */
         DynamoTraits::CodType   mCodC;

         /**
          * @brief Magnetic field
          */
         DynamoTraits::MagType   mMagB;

         /**
          * @brief Velocity field
          */
         DynamoTraits::VelType   mVelV;

         /**
          * @brief Induction equation
          */
         InductionMHD<DynamoTraits>    mInduction;

         /**
          * @brief Transport equation
          */
         TransportMHD<DynamoTraits>    mTransport;

         /**
          * @brief Navier Stokes equation
          */
         NavierStokesMHD<DynamoTraits>    mNavierStokes;
   };

}

#endif // DYNAMOSIMULATION_HPP
