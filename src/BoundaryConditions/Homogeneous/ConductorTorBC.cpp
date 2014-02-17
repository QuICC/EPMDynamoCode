/** \file ConductorTorBC.cpp
 *  \brief Implementation of the magnetic finite conductor boundary condition
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "BoundaryConditions/Homogeneous/ConductorTorBC.hpp"

// Project includes
//
#include "IO/ASCII/ConfigurationFile.hpp"

namespace EPMDynamo {

   ConductorTorBC::ConductorTorBC(const ConductorTorBC::BasisType &basis)
      : HomogeneousBC(basis.basisN(), basis.polyN())
   {
      // Setup the parametrisation
      this->setup();

      // Fill Operator BC values
      this->fillLHSBCValues(basis);
   }

   void ConductorTorBC::setup()
   {
      std::vector<std::string>   integers;
      std::vector<std::string>   floats;

      floats.push_back("cprime");

      ConfigurationFile cfg("conductor", integers, floats);

      // Initialise configuration file
      cfg.init();

      // Read data
      cfg.read();

      // Finalise configuration file reader
      cfg.finalise();

      // Print information
      cfg.printInfo();

      // Store values for configuration file
      this->mCprime = cfg.floats()(0);
   }

   void ConductorTorBC::fillLHSBCValues(const ConductorTorBC::BasisType &basis)
   {
      int nL = this->nL();
      for(int l = 0; l < nL; ++l)
      {
         this->rLHSBCValues(l) = (1 + this->mCprime)*basis.at(l).bpoly() + this->mCprime*basis.at(l).bdiff(1);
      }
   }

}
