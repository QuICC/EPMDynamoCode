/** \file SimulationConstants.hpp
 *  \brief Some useful simulation wide static constants
 */

#ifndef SIMULATIONCONSTANTS_HPP
#define SIMULATIONCONSTANTS_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "Polynomials/WorlandChebyshevTraits.hpp"
#include "Polynomials/WorlandLegendreTraits.hpp"

namespace EPMDynamo {

   /**
    * \brief Some useful simulation wide static constants
    */
   class SimulationConstants
   {
      public:
         /// Typedef for the radial worland polynomials traits to generate (either Chebyshev grid or Legendre grid)
           typedef WorlandChebyshevTraits RadialGridTraits;
         // typedef WorlandLegendreTraits RadialGridTraits;

         /**
          * @brief Check if grid is a CSCS grid
          */
         static bool isCSCSGrid();

         /**
          * @brief Test if grid type is locked
          */
         static bool isGridLocked();

         /**
          * @brief Activate CSCS grid computation
          */
         static void setCSCSGrid();

         /**
          * @brief Lock the grid type setting
          */
         static void lockGridType();

         /**
          * @brief Simple empty destructor
          */
         virtual ~SimulationConstants() {};

      protected:
         /**
          * @brief Flag for CSCS grid type
          */
         static bool msIsCSCS;

         /**
          * @brief Selection of grid type is locked?
          */
         static bool msIsGridLocked;

         /**
          * @brief Constructor
          */
         SimulationConstants() {};

      private:
   };

   inline bool SimulationConstants::isCSCSGrid()
   {
      return SimulationConstants::msIsCSCS;
   }

   inline bool SimulationConstants::isGridLocked()
   {
      return SimulationConstants::msIsGridLocked;
   }
}

#endif // SIMULATIONCONSTANTS_HPP
