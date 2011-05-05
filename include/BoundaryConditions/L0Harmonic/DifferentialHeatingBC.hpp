/** \file DifferentialHeatingBC.hpp
 *  \brief Implementation of the l=0 differential heating boundary condition
 */

#ifndef DIFFERENTIALHEATINGBC_HPP
#define DIFFERENTIALHEATINGBC_HPP

// Configuration includes
//
#include "Config/SimulationConfig.hpp"

// System includes
//

// External includes
//

// Project includes
//
#include "General/EPMTypedefs.hpp"
#include "BoundaryConditions/L0HarmonicBC.hpp"

namespace EPMDynamo {

   /**
    * @brief This class defines the  L=0 harmonic differential heating boundary condition
    */
   class DifferentialHeatingBC: public L0HarmonicBC
   {
      public:
         /// Typedef from Simulation trait to local truncation type
         typedef SimulationConfig::NumericalScheme::RadialBasisType    BasisType;

         /**
          * @brief Constructor
          *
          * @param bcValue Value at the boundary
          * @param basis   Radial basis
          */
         DifferentialHeatingBC(const EPMFloat bcValue, const BasisType &basis);

         /**
          * @brief Destructor
          */
         virtual ~DifferentialHeatingBC() {};

      protected:

         /**
          * @brief Fill in values for the boundaries
          *
          * @param basis   Radial basis
          */
         void fillLHSBCValues(const BasisType &basis);

         /**
          * @brief Set Boundary value
          *
          * @param bcValue   Value at the boundary
          */
         void setBCValue(const EPMFloat bcValue);

      private:
   };

}

#endif // DIFFERENTIALHEATINGBC_HPP
