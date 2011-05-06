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
         
      protected:

      private:
   };

   inline SmartTruncation WSHSimulation::createTrunc(const int maxN, const int maxL, const int maxM, const int mp, const int nCore)
   {
      // Lock grid type
      SimulationConstants::lockGridType();
      
      // Set radial grid size
      int nR = std::ceil(3.0*maxN/2.0) + std::ceil(3.0*maxL/4.0) + 1;

      return SmartTruncation(new Truncation(maxN, nR, maxL, maxM, mp, nCore));
   }

   inline SmartTruncation WSHSimulation::createTrunc(const ArrayI& arr)
   {
      // Lock grid type
      SimulationConstants::lockGridType();

      // Set radial grid size
      int nR = std::ceil(3.0*arr(0)/2.0) + std::ceil(3.0*arr(1)/4.0) + 1;

      return SmartTruncation(new Truncation(arr(0), nR, arr(1), arr(2), arr(3), arr(4)));
   }
}

#endif // WSHSIMULATION_HPP
