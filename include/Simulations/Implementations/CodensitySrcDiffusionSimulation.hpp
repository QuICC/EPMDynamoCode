/** \file CodensitySrcDiffusionSimulation.hpp
 *  \brief Implementation of a codensity diffusion simulation
 */

#ifndef CODENSITYSRCDIFFUSIONSIMULATION_HPP
#define CODENSITYSRCDIFFUSIONSIMULATION_HPP

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
#include "Simulations/Traits/CodensitySrcDiffusionTraits.hpp"

#include "IO/HDF5/State/StateFileReader.hpp"
#include "IO/HDF5/State/StateFileWriter.hpp"
#include "IO/HDF5/Source/SourceFileReader.hpp"
#include "IO/ASCII/EnergyFile.hpp"
#include "IO/ASCII/SpectrumFile.hpp"
#include "IO/ASCII/TimeFile.hpp"

#include "Equations/Transport/TransportSource.hpp"

#include "BoundaryConditions/Homogeneous/ZeroBC.hpp"

namespace EPMDynamo {

   /**
    * \brief Implementation of a codensity diffusion simulation
    */
   class CodensitySrcDiffusionSimulation: public SimulationBase
   {
      public:
         /**
          * @brief Simple empty destructor
          */
         virtual ~CodensitySrcDiffusionSimulation() {};

      protected:
         /**
          * @brief Constructor
          */
         CodensitySrcDiffusionSimulation();

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
         CodensitySrcDiffusionTraits::CodType   mCodC;

         /**
          * @brief Transport equation
          */
         TransportSource<CodensitySrcDiffusionTraits>    mTransport;
   };

}

#endif // CODENSITYSRCDIFFUSIONSIMULATION_HPP
