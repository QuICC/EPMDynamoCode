/** \file ConductorPolBC.cpp
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
#include "BoundaryConditions/Homogeneous/ConductorPolBC.hpp"

// Project includes
//
#include "IO/ASCII/ConfigurationFile.hpp"

namespace EPMDynamo {

   ConductorPolBC::ConductorPolBC(const ConductorPolBC::BasisType &basis)
      : HomogeneousBC(basis.basisN(), basis.polyN())
   {
      // Setup the parametrisation
      this->setup();

      // Fill Operator BC values
      this->fillLHSBCValues(basis);
   }

   void ConductorPolBC::setup()
   {
      std::vector<std::string>   integers;
      std::vector<std::string>   floats;

      floats.push_back("muratio");
      floats.push_back("c");
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
      this->mMuratio = cfg.floats()(0);
      this->mC = cfg.floats()(1);
      this->mCprime = cfg.floats()(2);
   }

   void ConductorPolBC::fillLHSBCValues(const ConductorPolBC::BasisType &basis)
   {
      EPMFloat rl = 0;

      int nL = this->nL();
      for(int l = 0; l < nL; ++l)
      {
         rl = static_cast<EPMFloat>(basis.at(l).l());

         this->rLHSBCValues(l) = (this->mMuratio + rl)*basis.at(l).bpoly() + this->mMuratio*basis.at(l).bdiff(1);

         this->rLHSBCValues(l) += (this->mC*rl)*(basis.at(l).bpoly() + basis.at(l).bdiff(1));

         this->rLHSBCValues(l) += this->mCprime*(this->mMuratio + this->mC*rl/2.0)*(2.0*basis.at(l).bdiff(1) + basis.at(l).bdiff(2) - rl*(rl+1.0)*basis.at(l).bpoly());
      }
   }

}
