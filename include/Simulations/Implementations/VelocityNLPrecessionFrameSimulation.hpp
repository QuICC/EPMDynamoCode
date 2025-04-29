/** \file VelocityNLPrecessionFrameSimulation.hpp
 *  \brief Implementation of a velocity precession frame simulation
 */

#ifndef VELOCITYNLPRECESSIONFRAMESIMULATION_HPP
#define VELOCITYNLPRECESSIONFRAMESIMULATION_HPP

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
#include "Simulations/Traits/VelocityNLRotDiffusionTraits.hpp"

#include "IO/HDF5/State/StateFileReader.hpp"
#include "IO/HDF5/State/StateFileWriter.hpp"
#include "IO/HDF5/Source/SourceFileReader.hpp"
#include "IO/ASCII/EnergyFile.hpp"
#include "IO/ASCII/SpectrumFile.hpp"
#include "IO/ASCII/TimeFile.hpp"

#include "Equations/NavierStokes/NavierStokesNLPrecessionFrame.hpp"

#include "BoundaryConditions/Homogeneous/ZeroBC.hpp"
#include "BoundaryConditions/Inhomogeneous/PrecessionFrameBC.hpp"
#include "BoundaryConditions/Homogeneous/DDRadialBC.hpp"
#include "BoundaryConditions/Homogeneous/DRadialBC.hpp"
#include "BoundaryConditions/Homogeneous/StressFreeTorBC.hpp"

namespace EPMDynamo {

   /**
    * \brief Implementation of a velocity diffusion simulation
    */
   class VelocityNLPrecessionFrameSimulation: public SimulationBase
   {
      public:
         /**
          * @brief Simple empty destructor
          */
         virtual ~VelocityNLPrecessionFrameSimulation() {};

      protected:
         /**
          * @brief Constructor
          */
         VelocityNLPrecessionFrameSimulation();

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
          * @brief Velocity field
          */
         VelocityNLRotDiffusionTraits::VelType   mVelV;

         /**
          * @brief Navier Stokes equation
          */
         NavierStokesNLPrecessionFrame<VelocityNLRotDiffusionTraits>    mNavierStokes;
   };

}

#endif // VELOCITYNLPRECESSIONFRAMESIMULATION_HPP
