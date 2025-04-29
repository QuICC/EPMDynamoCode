/** \file EEkEmRaParameters.cpp
 *  \brief Source of the implementation of the equation non dimensional parameters
 */

// System includes
//

// External includes
//

// Class include
//
#include "Equations/Parameters/EEkEmRaParameters.hpp"

// Project includes
//

namespace EPMDynamo {

   const std::string EEkEmRaParameters::id = "EEkEmRa";

   EEkEmRaParameters::EEkEmRaParameters(EPMFloat E, EPMFloat Ek, EPMFloat Em, EPMFloat Ra)
      : EquationParameters(), mE(E), mEk(Ek), mEm(Em), mRa(Ra)
   {
   }

   EEkEmRaParameters::EEkEmRaParameters(const Array& arr)
      : EquationParameters(), mE(arr(0)), mEk(arr(1)), mEm(arr(2)), mRa(arr(3))
   {
   }

   Array EEkEmRaParameters::params() const
   {
      Array values(this->nParams());

      values(0) = this->E();
      values(1) = this->Ek();
      values(2) = this->Em();
      values(3) = this->Ra();

      return values;
   }

   EPMFloat EEkEmRaParameters::nsDt() const
   {
      return 1.0;
   }

   EPMFloat EEkEmRaParameters::nsDiffusion() const
   {
      return this->E();
   }

   EPMFloat EEkEmRaParameters::nsLorentz() const
   {
      return 1.0;
   }

   EPMFloat EEkEmRaParameters::nsAdvection() const
   {
      return 1.0;
   }

   EPMFloat EEkEmRaParameters::nsBuoyancy() const
   {
      return this->Ra();
   }

   EPMFloat EEkEmRaParameters::nsCoriolis() const
   {
      return 2.0;
   }

   EPMFloat EEkEmRaParameters::indDt() const
   {
      return 1.0;
   }

   EPMFloat EEkEmRaParameters::indDiffusion() const
   {
      return this->Em();
   }

   EPMFloat EEkEmRaParameters::indAdvection() const
   {
      return 1.0;
   }

   EPMFloat EEkEmRaParameters::tptDt() const
   {
      return 1.0;
   }

   EPMFloat EEkEmRaParameters::tptDiffusion() const
   {
      return this->Ek();
   }

   EPMFloat EEkEmRaParameters::tptAdvection() const
   {
      return 1.0;
   }

   EPMFloat EEkEmRaParameters::alfvenFactor() const
   {
      return 1.0/(this->Em()*this->Em());
   }

   EPMFloat EEkEmRaParameters::alfvenDamping(const EPMFloat delta) const
   {
      return std::pow((this->E()+this->Em())/(2.0*delta*this->Em()),2);
   }

   EPMFloat EEkEmRaParameters::keFactor() const
   {
      return 1.0;
   }

   EPMFloat EEkEmRaParameters::meFactor() const
   {
      return 1.0;
   }

   EPMFloat EEkEmRaParameters::inertialCFL() const
   {
      return 1.0;
   }

   EPMFloat EEkEmRaParameters::torsionalCFL() const
   {
      return 1.0;
   }

   EPMFloat EEkEmRaParameters::codSourceScale() const
   {
      return 1.0;
   }

   EPMFloat EEkEmRaParameters::imposedMagScale() const
   {
      return 1.0;
   }

}
