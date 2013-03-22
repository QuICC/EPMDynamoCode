/** \file ONParameters.cpp
 *  \brief Source of the implementation of the equation non dimensional parameters
 */

// System includes
//

// External includes
//

// Class include
//
#include "Equations/Parameters/ONParameters.hpp"

// Project includes
//

namespace EPMDynamo {

   const std::string ONParameters::id = "ON";

   ONParameters::ONParameters(EPMFloat omega, EPMFloat nu)
      : EquationParameters(), mOmega(omega), mNu(nu)
   {
   }

   ONParameters::ONParameters(const Array& arr)
      : EquationParameters(), mOmega(arr(0)), mNu(arr(1))
   {
   }

   Array ONParameters::params() const
   {
      Array values(this->nParams());

      values(0) = this->omega();
      values(1) = this->nu();

      return values;
   }

   EPMFloat ONParameters::nsDt() const
   {
      return 1.0;
   }

   EPMFloat ONParameters::nsDiffusion() const
   {
      return this->nu();
   }

   EPMFloat ONParameters::nsLorentz() const
   {
      return 0.0;
   }

   EPMFloat ONParameters::nsAdvection() const
   {
      return 1.0;
   }

   EPMFloat ONParameters::nsBuoyancy() const
   {
      return 0.0;
   }

   EPMFloat ONParameters::nsCoriolis() const
   {
      return 2.0*this->omega();
   }

   EPMFloat ONParameters::indDt() const
   {
      return 0.0;
   }

   EPMFloat ONParameters::indDiffusion() const
   {
      return 0.0;
   }

   EPMFloat ONParameters::indAdvection() const
   {
      return 0.0;
   }

   EPMFloat ONParameters::tptDt() const
   {
      return 0.0;
   }

   EPMFloat ONParameters::tptDiffusion() const
   {
      return 0.0;
   }

   EPMFloat ONParameters::tptAdvection() const
   {
      return 0.0;
   }

   EPMFloat ONParameters::alfvenFactor() const
   {
      return 0.0;
   }

   EPMFloat ONParameters::alfvenDamping(const EPMFloat delta) const
   {
      return 0.0;
   }

   EPMFloat ONParameters::keFactor() const
   {
      return 0.5;
   }

   EPMFloat ONParameters::meFactor() const
   {
      return 0.0;
   }

   EPMFloat ONParameters::inertialCFL() const
   {
      // 1/2 * Ekman
      return (0.02)*this->nu()/(2.*this->omega());
   }

   EPMFloat ONParameters::torsionalCFL() const
   {
      return 1.0;
   }

   EPMFloat ONParameters::codSourceScale() const
   {
      return 0.0;
   }

   EPMFloat ONParameters::imposedMagScale() const
   {
      return 0.0;
   }

}
