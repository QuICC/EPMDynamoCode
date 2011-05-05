/** \file ZeroBC.hpp
 *  \brief Implementation of the basis Zero value boundary condition
 */

#ifndef ZEROBC_HPP
#define ZEROBC_HPP

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
    * @brief This class defines the zero boundary condition
    */
   class ZeroBC: public HomogeneousBC
   {
      public:
         /// Typedef from Simulation trait to local truncation type
         typedef SimulationConfig::NumericalScheme::RadialBasisType    BasisType;

         /**
          * @brief Constructor
          *
          * @param basis   Radial basis
          */
         ZeroBC(const BasisType &basis);

         /**
          * @brief Destructor
          */
         virtual ~ZeroBC() {};

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

#endif // ZEROBC_HPP
