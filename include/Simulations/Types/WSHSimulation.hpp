/** \file WSHSimulation.hpp
 *  \brief Definition of a simulation type using Worland radial expansion and spherical harmonics
 */

#ifndef WSHSIMULATION_HPP
#define WSHSIMULATION_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "GeneralScalars/SpectralSHScalar.hpp"
#include "Polynomials/RadialBasis.hpp"
#include "Polynomials/WorlandPolynomial.hpp"
#include "Transforms/Polynomial/RadialTransform.hpp"
#include "Simulations/SimulationConstants.hpp"
#include "GridGenerators/JacobiGrid.hpp"

namespace EPMDynamo {

   /**
    * \brief Definition of a simulation type using Worland radial expansion and spherical harmonics
    */
   class WSHSimulation
   {
      public:
         /// Typedef for the scalar data type
         typedef SpectralSHScalar ScalarType;

         /// Typedef for the radial basis data type
         typedef RadialBasis<WorlandPolynomial> RadialBasisType;

         /// Typedef for the radial transform data type
         typedef RadialTransform<WorlandPolynomial> RadialTransformType;

         /**
          * @brief Constructor
          */
         WSHSimulation();

         /**
          * @brief Destructor
          */
         virtual ~WSHSimulation() {};

         /**
          * @brief Create WSH simulation specific Spectral truncation object
          *
          * @param maxN Highest radial truncation
          * @param maxL Highest harmonic degree
          * @param maxM Highest harmonic order
          * @param mp  Longitudinal symmetry assumption
          * @param nCore Number of cores available
          */
         static SmartTruncation  createTrunc(const int maxN, const int maxL, const int maxM, const int mp, const int nCore);

         /**
          * @brief Create WSH simulation specific Spectral truncation object from parameters packed in array
          *
          * @param arr Array of truncation paramters (see above)
          */
         static SmartTruncation  createTrunc(const ArrayI &arr);

         /**
          * @brief Create radial spherical geometry factors
          */
         static SmartArray createRadSll(const int maxN, const int maxL, const Array& radGrid);
         
      protected:

         /**
          * @brief Compute the physical radial grid size
          *
          * @param maxN Highest radial truncation
          * @param maxL Highest harmonic degree
          */
         static int radialGridSize(const int maxN, const int maxL);

      private:
   };

   inline int WSHSimulation::radialGridSize(const int maxN, const int maxL)
   {
      // Set radial grid size. Aliasing requires to resolve (maxN+1)
      return std::ceil(3.0*(maxN + 1)/2.0) + std::ceil(3.0*maxL/4.0) + 1;
   }

   inline SmartTruncation WSHSimulation::createTrunc(const int maxN, const int maxL, const int maxM, const int mp, const int nCore)
   {
      // Lock grid type
      SimulationConstants::lockGridType();
      
      // Set radial grid size. Aliasing requires to resolve (maxN+1)
      int nR = WSHSimulation::radialGridSize(maxN, maxL);

      // Return the a smart truncation
      return SmartTruncation(new Truncation(maxN, nR, maxL, maxM, mp, nCore));
   }

   inline SmartTruncation WSHSimulation::createTrunc(const ArrayI& arr)
   {
      // Return the a smart truncation
      return WSHSimulation::createTrunc(arr(0), arr(1), arr(2), arr(3), arr(4));
   }

   inline SmartArray WSHSimulation::createRadSll(const int maxN, const int maxL, const Array& radGrid)
   {
      // Create the grid generator
      int gridN = static_cast<int>(std::ceil(static_cast<EPMFloat>(3*(maxN+1))/2.0));

      // Create the jacobi grid generator
      JacobiGrid   generator(gridN, GridGeneratorBase::RADIAL_GRID, SimulationConstants::isCSCSGrid());

      // Create storage for the jacobi grid
      Array jacobiGrid(maxN + 1);

      // Create the radSll smart array
      SmartArray  pRadSll = SmartArray(new Array(WSHSimulation::radialGridSize(maxN, maxL)));

      // Loop over all harmonic degrees starting from the highest
      for(int l = maxL; l > 0; --l)
      {
         // Recompute the grid for given parameters
         generator.recompute(SimulationConstants::RadialGridTraits::ALPHA, static_cast<EPMFloat>(l + SimulationConstants::RadialGridTraits::BETA));

         // Get the computed grid
         jacobiGrid = generator.grid();

         // Loop over the radiuses
         int r = 0;
         while(jacobiGrid(0) > radGrid(r))
         {
            r++;
            (*pRadSll)(r) = std::sqrt(static_cast<EPMFloat>(l*(l+1)));
         }

         // Special case for l=0
         (*pRadSll)(0) = std::sqrt(2.0);
      }

      return pRadSll;
   }
}

#endif // WSHSIMULATION_HPP
