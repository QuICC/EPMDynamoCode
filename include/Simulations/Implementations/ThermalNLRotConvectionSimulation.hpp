/** \file ThermalNLRotConvectionSimulation.hpp
 *  \brief Implementation of a thermal convection simulation
 */

#ifndef THERMALNLROTCONVECTIONSIMULATION_HPP
#define THERMALNLROTCONVECTIONSIMULATION_HPP

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
#include "Simulations/Traits/ThermalNLRotConvectionTraits.hpp"

#include "IO/HDF5/State/StateFileReader.hpp"
#include "IO/HDF5/State/StateFileWriter.hpp"
#include "IO/HDF5/Source/SourceFileReader.hpp"
#include "IO/ASCII/EnergyFile.hpp"
#include "IO/ASCII/FieldProbeFile.hpp"
#include "IO/ASCII/SpectrumFile.hpp"
#include "IO/ASCII/TimeFile.hpp"
#include "IO/ASCII/AngularMomentumFile.hpp"

#include "Equations/Transport/TransportMHD.hpp"
#include "Equations/NavierStokes/NavierStokesNLRotConvection.hpp"

#include "BoundaryConditions/Homogeneous/ZeroBC.hpp"
#include "BoundaryConditions/Homogeneous/DRadialBC.hpp"
#include "BoundaryConditions/Homogeneous/DDRadialBC.hpp"
#include "BoundaryConditions/Homogeneous/StressFreeTorBC.hpp"

namespace EPMDynamo {

   /**
    * \brief Implementation of a thermal convection simulation
    */
   class ThermalNLRotConvectionSimulation: public SimulationBase
   {
      public:
         /**
          * @brief Simple empty destructor
          */
         virtual ~ThermalNLRotConvectionSimulation() {};

      protected:
         /**
          * @brief Constructor
          */
         ThermalNLRotConvectionSimulation();

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
         ThermalNLRotConvectionTraits::CodType   mCodC;

         /**
          * @brief Velocity field
          */
         ThermalNLRotConvectionTraits::VelType   mVelV;

         /**
          * @brief Transport equation
          */
         TransportMHD<ThermalNLRotConvectionTraits>    mTransport;

         /**
          * @brief Navier Stokes equation
          */
         NavierStokesNLRotConvection<ThermalNLRotConvectionTraits>    mNavierStokes;
   };

}

#endif // THERMALNLROTCONVECTIONSIMULATION_HPP
