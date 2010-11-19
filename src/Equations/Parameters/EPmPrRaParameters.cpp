/** \file EPmPrRaParameters.cpp
 *  \brief Source of the implementation of the equation non dimensional parameters
 */

// System includes
//

// External includes
//

// Class include
//
#include "Equations/Parameters/EPmPrRaParameters.hpp"

// Project includes
//

namespace EPMDynamo {

   const std::string EPmPrRaParameters::id = "EPmPrRa";

   EPmPrRaParameters::EPmPrRaParameters(EPMFloat E, EPMFloat Pm, EPMFloat Pr, EPMFloat Ra)
      : EquationParameters(), mE(E), mPm(Pm), mPr(Pr), mRa(Ra)
   {
   }

   EPmPrRaParameters::EPmPrRaParameters(const Array& arr)
      : EquationParameters(), mE(arr(0)), mPm(arr(1)), mPr(arr(2)), mRa(arr(3))
   {
   }

   Array EPmPrRaParameters::params() const
   {
      Array values(this->nParams());

      values(0) = this->E();
      values(1) = this->Pm();
      values(2) = this->Pr();
      values(3) = this->Ra();

      return values;
   }

   EPMFloat EPmPrRaParameters::nsDt() const
   {
      return this->E();
   }

   EPMFloat EPmPrRaParameters::nsDiffusion() const
   {
      return this->E();
   }

   EPMFloat EPmPrRaParameters::nsLorentz() const
   {
      return 1.0/this->Pm();
   }

   EPMFloat EPmPrRaParameters::nsAdvection() const
   {
      return this->E();
   }

   EPMFloat EPmPrRaParameters::nsBuoyancy() const
   {
      return this->Ra();
   }

   EPMFloat EPmPrRaParameters::nsCoriolis() const
   {
      return 2.0;
   }

   EPMFloat EPmPrRaParameters::indDt() const
   {
      return 1.0;
   }

   EPMFloat EPmPrRaParameters::indDiffusion() const
   {
      return 1.0/this->Pm();
   }

   EPMFloat EPmPrRaParameters::indAdvection() const
   {
      return 1.0;
   }

   EPMFloat EPmPrRaParameters::tptDt() const
   {
      return 1.0;
   }

   EPMFloat EPmPrRaParameters::tptDiffusion() const
   {
      return 1.0/this->Pr();
   }

   EPMFloat EPmPrRaParameters::tptAdvection() const
   {
      return 1.0;
   }

   EPMFloat EPmPrRaParameters::alfvenFactor() const
   {
      return this->Pm()/this->E();
   }

   EPMFloat EPmPrRaParameters::alfvenDamping(const EPMFloat delta) const
   {
      return std::pow((this->Pm()+1)/(2.0*delta),2);
   }

   EPMFloat EPmPrRaParameters::meFactor() const
   {
      return 1.0;
   }

   void EPmPrRaParameters::testGlobalCFL(EPMFloat &rDt) const
   {
   }

   EPMFloat EPmPrRaParameters::codSourceScale() const
   {
      return 1.0;
   }

   EPMFloat EPmPrRaParameters::imposedMagScale() const
   {
      return 1.0;
   }

}
