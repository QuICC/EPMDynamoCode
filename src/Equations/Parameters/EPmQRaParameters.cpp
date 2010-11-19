/** \file EPmQRaParameters.cpp
 *  \brief Source of the implementation of the equation non dimensional parameters
 */

// System includes
//

// External includes
//

// Class include
//
#include "Equations/Parameters/EPmQRaParameters.hpp"

// Project includes
//

namespace EPMDynamo {

   const std::string EPmQRaParameters::id = "EPmQRa";

   EPmQRaParameters::EPmQRaParameters(EPMFloat E, EPMFloat Pm, EPMFloat q, EPMFloat Ra)
      : EquationParameters(), mE(E), mPm(Pm), mQ(q), mRa(Ra)
   {
   }

   EPmQRaParameters::EPmQRaParameters(const Array& arr)
      : EquationParameters(), mE(arr(0)), mPm(arr(1)), mQ(arr(2)), mRa(arr(3))
   {
   }

   Array EPmQRaParameters::params() const
   {
      Array values(this->nParams());

      values(0) = this->E();
      values(1) = this->Pm();
      values(2) = this->q();
      values(3) = this->Ra();

      return values;
   }

   EPMFloat EPmQRaParameters::nsDt() const
   {
      return 1.0/this->Pm();
   }

   EPMFloat EPmQRaParameters::nsDiffusion() const
   {
      return 1.0;
   }

   EPMFloat EPmQRaParameters::nsLorentz() const
   {
      return 1.0/this->E();
   }

   EPMFloat EPmQRaParameters::nsAdvection() const
   {
      return 1.0/this->Pm();
   }

   EPMFloat EPmQRaParameters::nsBuoyancy() const
   {
      return this->Ra();
   }

   EPMFloat EPmQRaParameters::nsCoriolis() const
   {
      return 1.0/this->E();
   }

   EPMFloat EPmQRaParameters::indDt() const
   {
      return 1.0;
   }

   EPMFloat EPmQRaParameters::indDiffusion() const
   {
      return 1.0;
   }

   EPMFloat EPmQRaParameters::indAdvection() const
   {
      return 1.0;
   }

   EPMFloat EPmQRaParameters::tptDt() const
   {
      return 1.0/this->q();
   }

   EPMFloat EPmQRaParameters::tptDiffusion() const
   {
      return 1.0;
   }

   EPMFloat EPmQRaParameters::tptAdvection() const
   {
      return 1.0/this->q();
   }

   EPMFloat EPmQRaParameters::alfvenFactor() const
   {
      return this->Pm()/this->E();
   }

   EPMFloat EPmQRaParameters::alfvenDamping(const EPMFloat delta) const
   {
      return std::pow((this->Pm()+1)/(2.0*delta),2);
   }

   EPMFloat EPmQRaParameters::meFactor() const
   {
      return 1.0;
   }

   void EPmQRaParameters::testGlobalCFL(EPMFloat &rDt) const
   {
   }

   EPMFloat EPmQRaParameters::codSourceScale() const
   {
      return 1.0;
   }

   EPMFloat EPmQRaParameters::imposedMagScale() const
   {
      return 1.0;
   }

}
