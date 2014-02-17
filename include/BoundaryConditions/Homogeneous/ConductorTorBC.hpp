/** \file ConductorTorBC.hpp
 *  \brief Implementation of the magnetic finite conductor toroidal boundary condition
 *  See Robert, Glatzmaier & Clune, 2010, GAFD, DOI: 10.1080/03091920903439753
 */

#ifndef CONDUCTORTORBC_HPP
#define CONDUCTORTORBC_HPP

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
   class ConductorTorBC: public HomogeneousBC
   {
      public:
         /// Typedef from Simulation trait to local truncation type
         typedef SimulationConfig::NumericalScheme::RadialBasisType    BasisType;

         /**
          * @brief Constructor
          *
          * @param basis   Radial basis
          */
         ConductorTorBC(const BasisType &basis);

         /**
          * @brief Destructor
          */
         virtual ~ConductorTorBC() {};

      protected:

         /**
          * @brief Setup the parametrisation
          */
         void setup();

         /**
          * @brief Fill in values for the boundaries
          *
          * @param basis   Radial basis
          */
         void fillLHSBCValues(const BasisType &basis);

      private:
         /**
          * @brief C prime parameter
          */
         EPMFloat mCprime;
   };

}

#endif // CONDUCTORTORBC_HPP
