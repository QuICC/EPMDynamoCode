/** \file PrRaParameters.cpp
 *  \brief Source of the implementation of the equation non dimensional parameters
 */

// System includes
//

// External includes
//

// Class include
//
#include "Equations/Parameters/PrRaParameters.hpp"

// Project includes
//

namespace EPMDynamo {

   const std::string PrRaParameters::id = "PrRa";

   PrRaParameters::PrRaParameters(EPMFloat Pr, EPMFloat Ra)
      : EquationParameters(), mPr(Pr), mRa(Ra)
   {
   }

   PrRaParameters::PrRaParameters(const Array& arr)
      : EquationParameters(), mPr(arr(0)), mRa(arr(1))
   {
   }

   Array PrRaParameters::params() const
   {
      Array values(this->nParams());

      values(0) = this->Pr();
      values(1) = this->Ra();

      return values;
   }

   EPMFloat PrRaParameters::nsDt() const
   {
      return 1.0;
   }

   EPMFloat PrRaParameters::nsDiffusion() const
   {
      return this->Pr();
   }

   EPMFloat PrRaParameters::nsLorentz() const
   {
      return 0.0;
   }

   EPMFloat PrRaParameters::nsAdvection() const
   {
      return 1.0;
   }

   EPMFloat PrRaParameters::nsBuoyancy() const
   {
      return this->Ra()/this->Pr();
   }

   EPMFloat PrRaParameters::nsCoriolis() const
   {
      return 0.0;
   }

   EPMFloat PrRaParameters::indDt() const
   {
      return 0.0;
   }

   EPMFloat PrRaParameters::indDiffusion() const
   {
      return 0.0;
   }

   EPMFloat PrRaParameters::indAdvection() const
   {
      return 0.0;
   }

   EPMFloat PrRaParameters::tptDt() const
   {
      return 1.0;
   }

   EPMFloat PrRaParameters::tptDiffusion() const
   {
      return 1.0;
   }

   EPMFloat PrRaParameters::tptAdvection() const
   {
      return 1.0;
   }

   EPMFloat PrRaParameters::alfvenFactor() const
   {
      return 0.0;
   }

   EPMFloat PrRaParameters::alfvenDamping(const EPMFloat delta) const
   {
      return 0.0;
   }

   EPMFloat PrRaParameters::keFactor() const
   {
      return 1.0;
   }

   EPMFloat PrRaParameters::meFactor() const
   {
      return 0.0;
   }

   EPMFloat PrRaParameters::inertialCFL() const
   {
      return 1.0;
   }

   EPMFloat PrRaParameters::torsionalCFL() const
   {
      return 1.0;
   }

   EPMFloat PrRaParameters::codSourceScale() const
   {
      return 1.0;
   }

   EPMFloat PrRaParameters::imposedMagScale() const
   {
      return 0.0;
   }

}
