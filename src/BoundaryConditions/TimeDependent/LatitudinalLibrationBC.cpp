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

      floats.push_back("epsilonx");
      floats.push_back("epsilony");
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
      this->mEpsilonx = cfg.floats()(0);
      this->mEpsilony = cfg.floats()(1);
      this->mFrequency = cfg.floats()(2);
   }

   EPMFloat LatitudinalLibrationBC::evolvingRealFactor() const
   {
      return this->mEpsilonx*std::cos(this->mFrequency*this->mrTSParams.time());
   }

   EPMFloat LatitudinalLibrationBC::evolvingImagFactor() const
   {
      //return 0.0;
      // By Yufeng Lin April 4 2016
      // Libration around the y axis u=epsy*sin(omgL*time)
      return this->mEpsilony*std::sin(this->mFrequency*this->mrTSParams.time());
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
      // by Yufeng Lin April 4 2016
      // libration around Y axis corresponds to imagenary part of T11
      if(l == 1 && m == 1)
      {
         return 0.5*this->evolvingImagFactor();
      } else
      {
      return 0.0;
      }
   }

   EPMComplex LatitudinalLibrationBC::getRHSBC(const int l, const int m) const
   {
      if(l == 1 && m == 1)
      {
         // return EPMComplex(-0.5*this->evolvingRealFactor(), 0.0);
         // by Yufeng Lin April 4 2016
         // Add imagenary part
         return EPMComplex(-0.5*this->evolvingRealFactor(),0.5*this->evolvingImagFactor());
      } else
      {
         return EPMComplex(0.0);
      }
   }
}
