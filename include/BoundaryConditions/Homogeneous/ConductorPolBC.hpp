/** \file ConductorPolBC.hpp
 *  \brief Implementation of the magnetic finite conductor poloidal boundary condition.
 *  See Robert, Glatzmaier & Clune, 2010, GAFD, DOI: 10.1080/03091920903439753
 */

#ifndef CONDUCTORPOLBC_HPP
#define CONDUCTORPOLBC_HPP

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
   class ConductorPolBC: public HomogeneousBC
   {
      public:
         /// Typedef from Simulation trait to local truncation type
         typedef SimulationConfig::NumericalScheme::RadialBasisType    BasisType;

         /**
          * @brief Constructor
          *
          * @param basis   Radial basis
          */
         ConductorPolBC(const BasisType &basis);

         /**
          * @brief Destructor
          */
         virtual ~ConductorPolBC() {};

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
          * @brief Mu ratio parameter
          */
         EPMFloat mMuratio;

         /**
          * @brief C parameter
          */
         EPMFloat mC;

         /**
          * @brief C prime parameter
          */
         EPMFloat mCprime;
   };

}

#endif // CONDUCTORPOLBC_HPP
