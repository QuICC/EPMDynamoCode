/** \file InvPrRaParameters.cpp
 *  \brief Source of the implementation of the equation non dimensional parameters
 */

// System includes
//

// External includes
//

// Class include
//
#include "Equations/Parameters/InvPrRaParameters.hpp"

// Project includes
//

namespace EPMDynamo {

   const std::string InvPrRaParameters::id = "PrRa";

   InvPrRaParameters::InvPrRaParameters(EPMFloat Pr, EPMFloat Ra)
      : EquationParameters(), mPr(Pr), mRa(Ra)
   {
   }

   InvPrRaParameters::InvPrRaParameters(const Array& arr)
      : EquationParameters(), mPr(arr(0)), mRa(arr(1))
   {
   }

   Array InvPrRaParameters::params() const
   {
      Array values(this->nParams());

      values(0) = this->Pr();
      values(1) = this->Ra();

      return values;
   }

   EPMFloat InvPrRaParameters::nsDt() const
   {
      return 1.0;
   }

   EPMFloat InvPrRaParameters::nsDiffusion() const
   {
      return 1.0;
   }

   EPMFloat InvPrRaParameters::nsLorentz() const
   {
      return 0.0;
   }

   EPMFloat InvPrRaParameters::nsAdvection() const
   {
      return 1.0;
   }

   EPMFloat InvPrRaParameters::nsBuoyancy() const
   {
      return this->Ra()*this->Pr();
   }

   EPMFloat InvPrRaParameters::nsCoriolis() const
   {
      return 0.0;
   }

   EPMFloat InvPrRaParameters::indDt() const
   {
      return 0.0;
   }

   EPMFloat InvPrRaParameters::indDiffusion() const
   {
      return 0.0;
   }

   EPMFloat InvPrRaParameters::indAdvection() const
   {
      return 0.0;
   }

   EPMFloat InvPrRaParameters::tptDt() const
   {
      return 1.0;
   }

   EPMFloat InvPrRaParameters::tptDiffusion() const
   {
      return 1.0/this->Pr();
   }

   EPMFloat InvPrRaParameters::tptAdvection() const
   {
      return 1.0;
   }

   EPMFloat InvPrRaParameters::alfvenFactor() const
   {
      return 0.0;
   }

   EPMFloat InvPrRaParameters::alfvenDamping(const EPMFloat delta) const
   {
      return 0.0;
   }

   EPMFloat InvPrRaParameters::keFactor() const
   {
      return 1.0;
   }

   EPMFloat InvPrRaParameters::meFactor() const
   {
      return 0.0;
   }

   EPMFloat InvPrRaParameters::inertialCFL() const
   {
      return 1.0;
   }

   EPMFloat InvPrRaParameters::torsionalCFL() const
   {
      return 1.0;
   }

   EPMFloat InvPrRaParameters::codSourceScale() const
   {
      return 1.0/this->Pr();
   }

   EPMFloat InvPrRaParameters::imposedMagScale() const
   {
      return 0.0;
   }

}
