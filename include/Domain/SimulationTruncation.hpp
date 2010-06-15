/** \file SimulationTruncation.hpp
 *  \brief Implementation of the global simulation truncation
 */

#ifndef SIMULATIONTRUNCATION_HPP
#define SIMULATIONTRUNCATION_HPP

// Configuration includes
//
#include "Config/SmartPointer.h"

// System includes
//

// External includes
//

// Project includes
//
#include "Domain/Simulation/RadialDomain.hpp"
#include "Domain/Simulation/HorizontalDomain.hpp"

namespace EPMDynamo {

   /**
    * \brief Implementation of the global simulation wide truncation information
    */
   class SimulationTruncation
   {
      public:
         /**
          * @brief Constructor
          *
          * \param maxN Highest radial basis truncation
          * \param nR   Number of real space radial grid points
          * \param maxL Highest harmonic degree
          * \param maxM Highest harmonic order
          * \param mp   Longitudinal symmetry assumption
          */
         SimulationTruncation(const int maxN, const int nR, const int maxL, const int maxM, const int mp);

         /**
          * @brief Destructor
          */
         virtual ~SimulationTruncation() {};

         /**
          * @brief Get the radial truncation information
          */
         const SmartSimRadDomain rad();

         /**
          * @brief Get the horizontal truncation information
          */
         const SmartSimHozDomain hoz();
         
      protected:

      private:

         /**
          * @brief Pointer to the radial truncation
          */
         SmartSimRadDomain mpRadial;

         /**
          * @brief Pointer to the horizontal truncation
          */
         SmartSimHozDomain mpHorizon;
   };

   inline const SmartSimRadDomain  SimulationTruncation::rad()
   {
      return this->mpRadial;
   }

   inline const SmartSimHozDomain  SimulationTruncation::hoz()
   {
      return this->mpHorizon;
   }

   /// Typedef for a smart SimulationTruncation
   typedef EPMSHARED_PTR<SimulationTruncation> SmartSimTrunc;

}

#endif // SIMULATIONTRUNCATION_HPP
