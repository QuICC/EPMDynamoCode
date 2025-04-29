/** \file InsulatingBC.hpp
 *  \brief Implementation of the magnetic insulating boundary condition
 */

#ifndef INSULATINGBC_HPP
#define INSULATINGBC_HPP

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
#include "BoundaryConditions/HomogeneousBC.hpp"

namespace EPMDynamo {

   /**
    * @brief This class defines the insulating homogeneous boundary condition
    */
   class InsulatingBC: public HomogeneousBC
   {
      public:
         /// Typedef from Simulation trait to local truncation type
         typedef SimulationConfig::NumericalScheme::RadialBasisType    BasisType;

         /**
          * @brief Constructor
          *
          * @param basis   Radial basis
          */
         InsulatingBC(const BasisType &basis);

         /**
          * @brief Destructor
          */
         virtual ~InsulatingBC() {};

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

#endif // INSULATINGBC_HPP
