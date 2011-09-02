/** \file LatitudinalLibrationBC.cpp
 *  \brief Implementation of the time dependent poloidal libration boundary condition
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "BoundaryConditions/TimeDependent/LatitudinalLibrationBC.hpp"

// Project includes
//
#include "IO/ASCII/ConfigurationFile.hpp"

namespace EPMDynamo {

   LatitudinalLibrationBC::LatitudinalLibrationBC(const LatitudinalLibrationBC::BasisType &basis, const TimestepParameters &tsParams)
      : TimeDependentBC(basis.basisN(), basis.polyN(), tsParams)
   {
      // Setup the parametrisation
      this->setup();

      // Fill Operator BC values
      this->fillLHSBCValues(basis);
   }

   void LatitudinalLibrationBC::setup()
   {
      std::vector<std::string>   integers;
      std::vector<std::string>   floats;

      floats.push_back("epsilon");
      floats.push_back("frequency");

      ConfigurationFile cfg("libration", integers, floats);

      // Initialise configuration file
      cfg.init();

      // Read data
      cfg.read();

      // Finalise configuration file reader
      cfg.finalise();

      // Print information
      cfg.printInfo();

      // Store values for configuration file
      this->mEpsilon = cfg.floats()(0);
      this->mFrequency = cfg.floats()(1);
   }

   EPMFloat LatitudinalLibrationBC::evolvingRealFactor() const
   {
      return this->mEpsilon*std::cos(this->mFrequency*this->mrTSParams.time());
   }

   EPMFloat LatitudinalLibrationBC::evolvingImagFactor() const
   {
      return 0.0;
   }

   void LatitudinalLibrationBC::fillLHSBCValues(const LatitudinalLibrationBC::BasisType &basis)
   {
      int nL = this->nL();
      for(int l = 0; l < nL; ++l)
      {
         rLHSBCValues(l) = basis.at(l).bpoly();
      }
   }

   EPMFloat LatitudinalLibrationBC::getRealRHSBC(const int l, const int m) const
   {
      if(l == 1 && m == 1)
      {
         return -0.5*this->evolvingRealFactor();
      } else
      {
         return 0.0;
      }
   }

   EPMFloat LatitudinalLibrationBC::getImagRHSBC(const int l, const int m) const
   {
      return 0.0;
   }

   EPMComplex LatitudinalLibrationBC::getRHSBC(const int l, const int m) const
   {
      if(l == 1 && m == 1)
      {
         return EPMComplex(-0.5*this->evolvingRealFactor(), 0.0);
      } else
      {
         return EPMComplex(0.0);
      }
   }
}
