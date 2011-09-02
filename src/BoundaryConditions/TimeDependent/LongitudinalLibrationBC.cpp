/** \file LongitudinalLibrationBC.cpp
 *  \brief Implementation of the time dependent toroidal longitudinal libration boundary condition
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "BoundaryConditions/TimeDependent/LongitudinalLibrationBC.hpp"

// Project includes
//
#include "IO/ASCII/ConfigurationFile.hpp"

namespace EPMDynamo {

   LongitudinalLibrationBC::LongitudinalLibrationBC(const LongitudinalLibrationBC::BasisType &basis, const TimestepParameters &tsParams)
      : TimeDependentBC(basis.basisN(), basis.polyN(), tsParams)
   {
      // Setup the parametrisation
      this->setup();

      // Fill Operator BC values
      this->fillLHSBCValues(basis);
   }

   void LongitudinalLibrationBC::setup()
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

   EPMFloat LongitudinalLibrationBC::evolvingRealFactor() const
   {
      return this->mEpsilon*std::cos(this->mFrequency*this->mrTSParams.time());
   }

   EPMFloat LongitudinalLibrationBC::evolvingImagFactor() const
   {
      return 0.0;
   }

   void LongitudinalLibrationBC::fillLHSBCValues(const LongitudinalLibrationBC::BasisType &basis)
   {
      int nL = this->nL();
      for(int l = 0; l < nL; ++l)
      {
         rLHSBCValues(l) = basis.at(l).bpoly();
      }
   }

   EPMFloat LongitudinalLibrationBC::getRealRHSBC(const int l, const int m) const
   {
      if(l == 1 && m == 0)
      {
         return -1.0*this->evolvingRealFactor();
      } else
      {
         return 0.0;
      }
   }

   EPMFloat LongitudinalLibrationBC::getImagRHSBC(const int l, const int m) const
   {
      return 0.0;
   }

   EPMComplex LongitudinalLibrationBC::getRHSBC(const int l, const int m) const
   {
      if(l == 1 && m == 0)
      {
         return EPMComplex(-1.0*this->evolvingRealFactor(), 0.0);
      } else
      {
         return EPMComplex(0.0);
      }
   }
}
