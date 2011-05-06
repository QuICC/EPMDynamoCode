/** \file StressFreeTorBC.hpp
 *  \brief Implementation of the stress free velocity boundary condition on the toroidal scalar
 */

#ifndef STRESSFREETORBC_HPP
#define STRESSFREETORBC_HPP

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
   class StressFreeTorBC: public HomogeneousBC
   {
      public:
         /// Typedef from Simulation trait to local truncation type
         typedef SimulationConfig::NumericalScheme::RadialBasisType    BasisType;

         /**
          * @brief Constructor
          *
          * @param basis   Radial basis
          */
         StressFreeTorBC(const BasisType &basis);

         /**
          * @brief Destructor
          */
         virtual ~StressFreeTorBC() {};

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

#endif // STRESSFREETORBC_HPP
