/** \file ELPmPrRaParameters.cpp
 *  \brief Source of the implementation of the equation non dimensional parameters
 */

// System includes
//

// External includes
//

// Class include
//
#include "Equations/Parameters/ELPmPrRaParameters.hpp"

// Project includes
//

namespace EPMDynamo {

   const std::string ELPmPrRaParameters::id = "ELPmPrRa";

   ELPmPrRaParameters::ELPmPrRaParameters(EPMFloat E, EPMFloat L, EPMFloat Pm, EPMFloat Pr, EPMFloat Ra)
      : EquationParameters(), mE(E), mL(L), mPm(Pm), mPr(Pr), mRa(Ra)
   {
   }

   ELPmPrRaParameters::ELPmPrRaParameters(const Array& arr)
      : EquationParameters(), mE(arr(0)), mL(arr(1)), mPm(arr(2)), mPr(arr(3)), mRa(arr(4))
   {
   }

   Array ELPmPrRaParameters::params() const
   {
      Array values(this->nParams());

      values(0) = this->E();
      values(1) = this->L();
      values(2) = this->Pm();
      values(3) = this->Pr();
      values(4) = this->Ra();

      return values;
   }

   EPMFloat ELPmPrRaParameters::nsDt() const
   {
      return this->E();
   }

   EPMFloat ELPmPrRaParameters::nsDiffusion() const
   {
      return this->E();
   }

   EPMFloat ELPmPrRaParameters::nsLorentz() const
   {
      return this->L()*this->Pm();
   }

   EPMFloat ELPmPrRaParameters::nsAdvection() const
   {
      return this->E();
   }

   EPMFloat ELPmPrRaParameters::nsBuoyancy() const
   {
      return this->E()*this->Ra();
   }

   EPMFloat ELPmPrRaParameters::nsCoriolis() const
   {
      return 1.0;
   }

   EPMFloat ELPmPrRaParameters::indDt() const
   {
      return 1.0;
   }

   EPMFloat ELPmPrRaParameters::indDiffusion() const
   {
      return 1.0/this->Pm();
   }

   EPMFloat ELPmPrRaParameters::indAdvection() const
   {
      return 1.0;
   }

   EPMFloat ELPmPrRaParameters::tptDt() const
   {
      return 1.0;
   }

   EPMFloat ELPmPrRaParameters::tptDiffusion() const
   {
      return 1.0/this->Pr();
   }

   EPMFloat ELPmPrRaParameters::tptAdvection() const
   {
      return 1.0;
   }

   EPMFloat ELPmPrRaParameters::alfvenFactor() const
   {
      return this->L()*this->Pm()/this->E();
   }

   EPMFloat ELPmPrRaParameters::alfvenDamping(const EPMFloat delta) const
   {
      return std::pow((1.0+1.0/this->Pm())/(2.0*delta),2);
   }

   EPMFloat ELPmPrRaParameters::keFactor() const
   {
      return this->Pm()*this->Pm();
   }

   EPMFloat ELPmPrRaParameters::meFactor() const
   {
      return this->L()*std::pow(this->Pm(),3)/this->E();
   }

   void ELPmPrRaParameters::testGlobalCFL(EPMFloat &rDt) const
   {
   }

   EPMFloat ELPmPrRaParameters::codSourceScale() const
   {
      return 1.0;
   }

   EPMFloat ELPmPrRaParameters::imposedMagScale() const
   {
      return 1.0/this->Pm();
   }

}
