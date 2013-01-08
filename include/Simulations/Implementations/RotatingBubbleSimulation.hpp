/** \file RotatingBubbleSimulation.hpp
 *  \brief Implementation of the rotating bubble simulation
 */

#ifndef ROTATINGBUBBLESIMULATION_HPP
#define ROTATINGBUBBLESIMULATION_HPP

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
#include "IO/ASCII/LibrationFile.hpp"
#include "IO/ASCII/AngularMomentumFile.hpp"

#include "Equations/NavierStokes/NavierStokesNLRotDiffusion.hpp"

#include "BoundaryConditions/Homogeneous/ZeroBC.hpp"
#include "BoundaryConditions/Homogeneous/DDRadialBC.hpp"
#include "BoundaryConditions/Homogeneous/DRadialBC.hpp"
#include "BoundaryConditions/Homogeneous/StressFreeTorBC.hpp"
#include "BoundaryConditions/Inhomogeneous/PrecessionFrameBC.hpp"
#include "BoundaryConditions/TimeDependent/LongitudinalLibrationBC.hpp"
#include "BoundaryConditions/TimeDependent/LatitudinalLibrationBC.hpp"
#include "BoundaryConditions/Inhomogeneous/RotatingBubbleNSBC.hpp"
#include "BoundaryConditions/Inhomogeneous/RotatingBubbleDRBC.hpp"

namespace EPMDynamo {

   /**
    * \brief Implementation of a rotating bubble simulation
    */
   class RotatingBubbleSimulation: public SimulationBase
   {
      public:
         /**
          * @brief Simple empty destructor
          */
         virtual ~RotatingBubbleSimulation() {};

      protected:
         /**
          * @brief Constructor
          */
         RotatingBubbleSimulation();

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
         NavierStokesNLRotDiffusion<VelocityNLRotDiffusionTraits>    mNavierStokes;
   };

}

#endif // ROTATINGBUBBLESIMULATION_HPP
