/** \file EParameters.cpp
 *  \brief Source of the implementation of the equation non dimensional parameters
 */

// System includes
//

// External includes
//

// Class include
//
#include "Equations/Parameters/EParameters.hpp"

// Project includes
//

namespace EPMDynamo {

   const std::string EParameters::id = "E";

   EParameters::EParameters(EPMFloat E)
      : EquationParameters(), mE(E)
   {
   }

   EParameters::EParameters(const Array& arr)
      : EquationParameters(), mE(arr(0))
   {
   }

   Array EParameters::params() const
   {
      Array values(this->nParams());

      values(0) = this->E();

      return values;
   }

   EPMFloat EParameters::nsDt() const
   {
      return 1.0;
   }

   EPMFloat EParameters::nsDiffusion() const
   {
      return this->E();
   }

   EPMFloat EParameters::nsLorentz() const
   {
      return 0.0;
   }

   EPMFloat EParameters::nsAdvection() const
   {
      return 1.0;
   }

   EPMFloat EParameters::nsBuoyancy() const
   {
      return 0.0;
   }

   EPMFloat EParameters::nsCoriolis() const
   {
      return 2.0;
   }

   EPMFloat EParameters::indDt() const
   {
      return 0.0;
   }

   EPMFloat EParameters::indDiffusion() const
   {
      return 0.0;
   }

   EPMFloat EParameters::indAdvection() const
   {
      return 0.0;
   }

   EPMFloat EParameters::tptDt() const
   {
      return 0.0;
   }

   EPMFloat EParameters::tptDiffusion() const
   {
      return 0.0;
   }

   EPMFloat EParameters::tptAdvection() const
   {
      return 0.0;
   }

   EPMFloat EParameters::alfvenFactor() const
   {
      return 0.0;
   }

   EPMFloat EParameters::alfvenDamping(const EPMFloat delta) const
   {
      return 0.0;
   }

   EPMFloat EParameters::keFactor() const
   {
      return 1.0;
   }

   EPMFloat EParameters::meFactor() const
   {
      return 0.0;
   }

   EPMFloat EParameters::inertialCFL() const
   {
      return 0.5/10.;
   }

   EPMFloat EParameters::torsionalCFL() const
   {
      return 1.0;
   }

   EPMFloat EParameters::codSourceScale() const
   {
      return 0.0;
   }

   EPMFloat EParameters::imposedMagScale() const
   {
      return 0.0;
   }

}
