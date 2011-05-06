/** \file MagnetoConvectionSimulation.hpp
 *  \brief Implementation of the magneto-convection simulation
 */

#ifndef MAGNETOCONVECTIONSIMULATION_HPP
#define MAGNETOCONVECTIONSIMULATION_HPP

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
#include "Simulations/Traits/MagnetoConvectionTraits.hpp"

#include "IO/HDF5/Imposed/ImposedFieldReader.hpp"
#include "IO/HDF5/State/StateFileReader.hpp"
#include "IO/HDF5/State/StateFileWriter.hpp"
#include "IO/HDF5/Source/SourceFileReader.hpp"
#include "IO/ASCII/EnergyFile.hpp"
#include "IO/ASCII/SpectrumFile.hpp"

#include "Equations/Induction/InductionMHD.hpp"
#include "Equations/NavierStokes/NavierStokesMHD.hpp"
#include "Equations/Transport/TransportMHD.hpp"

#include "BoundaryConditions/Homogeneous/ZeroBC.hpp"
#include "BoundaryConditions/Homogeneous/DDRadialBC.hpp"
#include "BoundaryConditions/Homogeneous/DRadialBC.hpp"
#include "BoundaryConditions/Homogeneous/InsulatingBC.hpp"
#include "BoundaryConditions/Homogeneous/StressFreeTorBC.hpp"

namespace EPMDynamo {

   /**
    * \brief Implementation of the magneto-convection simulation
    */
   class MagnetoConvectionSimulation: public SimulationBase
   {
      public:
         /**
          * @brief Simple empty destructor
          */
         virtual ~MagnetoConvectionSimulation() {};

      protected:
         /**
          * @brief Constructor
          */
         MagnetoConvectionSimulation();

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
         MagnetoConvectionTraits::CodType   mCodC;

         /**
          * @brief Magnetic field
          */
         MagnetoConvectionTraits::MagType   mMagB;

         /**
          * @brief Velocity field
          */
         MagnetoConvectionTraits::VelType   mVelV;

         /**
          * @brief Induction equation
          */
         InductionMHD<MagnetoConvectionTraits>    mInduction;

         /**
          * @brief Transport equation
          */
         TransportMHD<MagnetoConvectionTraits>    mTransport;

         /**
          * @brief Navier Stokes equation
          */
         NavierStokesMHD<MagnetoConvectionTraits>    mNavierStokes;
   };

}

#endif // MAGNETOCONVECTIONSIMULATION_HPP
