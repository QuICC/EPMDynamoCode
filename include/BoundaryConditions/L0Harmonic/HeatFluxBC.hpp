/** \file HeatFluxBC.hpp
 *  \brief Implementation of the l=0 constant heat flux boundary condition
 */

#ifndef HEATFLUXBC_HPP
#define HEATFLUXBC_HPP

// Configuration includes
//
#include "Config/SimulationConfig.hpp"
#include "Config/NumericalSchemeInc.hpp"

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
    * \brief Implementation of the l=0 constant heat flux boundary condition
    */
   class HeatFluxBC: public L0HarmonicBC
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
         HeatFluxBC(const EPMFloat bcValue, const BasisType &basis);

         /**
          * @brief Destructor
          */
         virtual ~HeatFluxBC() {};

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

#endif // HEATFLUXBC_HPP
