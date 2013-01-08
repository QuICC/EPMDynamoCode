/** \file StressFreeTorAngMomBC.hpp
 *  \brief Implementation of the stress free velocity boundary condition on the toroidal scalar including angular momentum conservation on l=1 modes.
 */

#ifndef STRESSFREETORANGMOMBC_HPP
#define STRESSFREETORANGMOMBC_HPP

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
    * @brief This class defines the stress free homogeneous boundary condition on a
    *    toroidal scalar
    */
   class StressFreeTorAngMomBC: public HomogeneousBC
   {
      public:
         /// Typedef from Simulation trait to local truncation type
         typedef SimulationConfig::NumericalScheme::RadialBasisType    BasisType;

         /**
          * @brief Constructor
          *
          * @param basis   Radial basis
          */
         StressFreeTorAngMomBC(const BasisType &basis);

         /**
          * @brief Destructor
          */
         virtual ~StressFreeTorAngMomBC() {};

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

#endif // STRESSFREETORANGMOMBC_HPP
