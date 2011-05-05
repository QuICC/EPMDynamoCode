/** \file DDRadialBC.hpp
 *  \brief Implementation of the homogeneous second radial derivative boundary condition
 */

#ifndef DDRADIALBC_HPP
#define DDRADIALBC_HPP

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
    * @brief This class defines the homogeneous boundary condition on the second radial
    *    derivative
    */
   class DDRadialBC: public HomogeneousBC
   {
      public:
         /// Typedef from Simulation trait to local truncation type
         typedef SimulationConfig::NumericalScheme::RadialBasisType    BasisType;

         /**
          * @brief Constructor
          *
          * @param basis   Radial basis
          */
         DDRadialBC(const BasisType &basis);

         /**
          * @brief Destructor
          */
         virtual ~DDRadialBC() {};

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

#endif // DDRADIALBC_HPP
