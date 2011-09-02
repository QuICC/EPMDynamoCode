/** \file EPmParameters.cpp
 *  \brief Source of the implementation of the equation non dimensional parameters
 */

// System includes
//

// External includes
//

// Class include
//
#include "Equations/Parameters/EPmParameters.hpp"

// Project includes
//

namespace EPMDynamo {

   const std::string EPmParameters::id = "EPm";

   EPmParameters::EPmParameters(EPMFloat E, EPMFloat Pm)
      : EquationParameters(), mE(E), mPm(Pm)
   {
   }

   EPmParameters::EPmParameters(const Array& arr)
      : EquationParameters(), mE(arr(0)), mPm(arr(1))
   {
   }

   Array EPmParameters::params() const
   {
      Array values(this->nParams());

      values(0) = this->E();
      values(1) = this->Pm();

      return values;
   }

   EPMFloat EPmParameters::nsDt() const
   {
      return 1.0;
   }

   EPMFloat EPmParameters::nsDiffusion() const
   {
      return this->E();
   }

   EPMFloat EPmParameters::nsLorentz() const
   {
      return 1.0;
   }

   EPMFloat EPmParameters::nsAdvection() const
   {
      return 1.0;
   }

   EPMFloat EPmParameters::nsBuoyancy() const
   {
      return 0.0;
   }

   EPMFloat EPmParameters::nsCoriolis() const
   {
      return 2.0;
   }

   EPMFloat EPmParameters::indDt() const
   {
      return 1.0;
   }

   EPMFloat EPmParameters::indDiffusion() const
   {
      return this->E()/this->Pm();
   }

   EPMFloat EPmParameters::indAdvection() const
   {
      return 1.0;
   }

   EPMFloat EPmParameters::tptDt() const
   {
      return 0.0;
   }

   EPMFloat EPmParameters::tptDiffusion() const
   {
      return 0.0;
   }

   EPMFloat EPmParameters::tptAdvection() const
   {
      return 0.0;
   }

   EPMFloat EPmParameters::alfvenFactor() const
   {
      return this->Pm()/this->E();
   }

   EPMFloat EPmParameters::alfvenDamping(const EPMFloat delta) const
   {
      return std::pow((this->Pm()+1)/(2.0*delta),2);
   }

   EPMFloat EPmParameters::keFactor() const
   {
      return 0.5;
   }

   EPMFloat EPmParameters::meFactor() const
   {
      return this->Pm()/this->E();
   }

   EPMFloat EPmParameters::inertialCFL() const
   {
      return 0.5/10.;
   }

   EPMFloat EPmParameters::torsionalCFL() const
   {
      return 1.0;
   }

   EPMFloat EPmParameters::codSourceScale() const
   {
      return 0.0;
   }

   EPMFloat EPmParameters::imposedMagScale() const
   {
      return 1.0;
   }

}
