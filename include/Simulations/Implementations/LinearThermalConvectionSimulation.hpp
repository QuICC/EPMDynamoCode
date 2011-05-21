/** \file LinearThermalConvectionSimulation.hpp
 *  \brief Implementation of a linear thermal convection simulation
 */

#ifndef LINEARTHERMALCONVECTIONSIMULATION_HPP
#define LINEARTHERMALCONVECTIONSIMULATION_HPP

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
#include "Simulations/Traits/LinearThermalConvectionTraits.hpp"

#include "IO/HDF5/State/StateFileReader.hpp"
#include "IO/HDF5/State/StateFileWriter.hpp"
#include "IO/HDF5/Source/SourceFileReader.hpp"
#include "IO/ASCII/EnergyFile.hpp"
#include "IO/ASCII/SpectrumFile.hpp"
#include "IO/ASCII/TimeFile.hpp"

#include "Equations/Transport/TransportLinear.hpp"
#include "Equations/NavierStokes/NavierStokesConvection.hpp"

#include "BoundaryConditions/Homogeneous/ZeroBC.hpp"
#include "BoundaryConditions/Homogeneous/DRadialBC.hpp"
#include "BoundaryConditions/Homogeneous/DDRadialBC.hpp"
#include "BoundaryConditions/Homogeneous/StressFreeTorBC.hpp"

namespace EPMDynamo {

   /**
    * \brief Implementation of a thermal convection simulation
    */
   class LinearThermalConvectionSimulation: public SimulationBase
   {
      public:
         /**
          * @brief Simple empty destructor
          */
         virtual ~LinearThermalConvectionSimulation() {};

      protected:
         /**
          * @brief Constructor
          */
         LinearThermalConvectionSimulation();

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
         LinearThermalConvectionTraits::CodType   mCodC;

         /**
          * @brief Velocity field
          */
         LinearThermalConvectionTraits::VelType   mVelV;

         /**
          * @brief Transport equation
          */
         TransportLinear<LinearThermalConvectionTraits>    mTransport;

         /**
          * @brief Navier Stokes equation
          */
         NavierStokesConvection<LinearThermalConvectionTraits>    mNavierStokes;
   };

}

#endif // LINEARTHERMALCONVECTIONSIMULATION_HPP
