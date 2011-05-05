/** \file DRadialBC.hpp
 *  \brief Implementation of the homogeneous first radial derivative boundary condition
 */

#ifndef DRADIALBC_HPP
#define DRADIALBC_HPP

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
#include "BoundaryConditions/HomogeneousBC.hpp"

namespace EPMDynamo {

   /**
    * @brief This class defines the homogeneous boundary condition on the first radial
    *    derivative
    */
   class DRadialBC: public HomogeneousBC
   {
      public:
         /// Typedef from Simulation trait to local truncation type
         typedef SimulationConfig::NumericalScheme::RadialBasisType    BasisType;

         /**
          * @brief Constructor
          *
          * @param basis   Radial basis
          */
         DRadialBC(const BasisType &basis);

         /**
          * @brief Destructor
          */
         virtual ~DRadialBC() {};

      protected:

         /**
          * @brief Fill in values for the boundaries
          *
          * @param basis   Radial basis
          */
         void fillLHSBCValues(const BasisType &basis);

      private:
   };
}

#endif // DRADIALBC_HPP
