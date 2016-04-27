/** \file ThermalPrecessionSimulation.hpp
 *  \brief Implementation of a thermal precession simulation
 */

#ifndef THERMALPRECESSIONFRAMESIMULATION_HPP
#define THERMALPRECESSIONFRAMESIMULATION_HPP

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
//#include "Simulations/Traits/ThermalPrecessionTraits.hpp"
#include "Simulations/Traits/ThermalNLRotConvectionTraits.hpp"

#include "IO/HDF5/State/StateFileReader.hpp"
#include "IO/HDF5/State/StateFileWriter.hpp"
#include "IO/HDF5/Source/SourceFileReader.hpp"
#include "IO/ASCII/EnergyFile.hpp"
#include "IO/ASCII/SpectrumFile.hpp"
#include "IO/ASCII/TimeFile.hpp"

//Leo: keeping this for induction 
//#include "Equations/Transport/InductionMHD.hpp"

#include "Equations/Transport/TransportMHD.hpp"
//#include "Equations/NavierStokes/NavierStokesRotConvection.hpp"
#include "Equations/NavierStokes/NavierStokesThermalPrecessionFrame.hpp"


#include "BoundaryConditions/Homogeneous/ZeroBC.hpp"
#include "BoundaryConditions/Homogeneous/DRadialBC.hpp"
#include "BoundaryConditions/Homogeneous/DDRadialBC.hpp"
#include "BoundaryConditions/Homogeneous/StressFreeTorBC.hpp"

namespace EPMDynamo {

   /**
    * \brief Implementation of a thermal precession frame simulation
    */
   class ThermalPrecessionFrameSimulation: public SimulationBase
   {
      public:
         /**
          * @brief Simple empty destructor
          */
         virtual ~ThermalPrecessionFrameSimulation() {};

      protected:
         /**
          * @brief Constructor
          */
         ThermalPrecessionFrameSimulation();

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
         //PrecessionDynamoTraits::MagType   mMagB;

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
         //InductionMHD<PrecessionDynamoTraits>    mInduction;


         /**
          * @brief Transport equation
          */
         TransportMHD<ThermalNLRotConvectionTraits>    mTransport;

         /**
          * @brief Navier Stokes equation
          */
         NavierStokesThermalPrecessionFrame<ThermalNLRotConvectionTraits> mNavierStokes;
   };

}

#endif // THERMALPRECESSIONFRAMESIMULATION_HPP
