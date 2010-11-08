/** \file EkQRaRoParameters.cpp
 *  \brief Source of the implementation of the equation non dimensional parameters
 */

// System includes
//

// External includes
//

// Class include
//
#include "Equations/Parameters/EkQRaRoParameters.hpp"

// Project includes
//

namespace EPMDynamo {

   EkQRaRoParameters::EkQRaRoParameters(EPMFloat E, EPMFloat q, EPMFloat Ra, EPMFloat Ro)
      : EquationParameters(), mE(E), mQ(q), mRa(Ra), mRo(Ro)
   {
   }

   EkQRaRoParameters::EkQRaRoParameters(const Array& arr)
      : EquationParameters(), mE(arr(0)), mQ(arr(1)), mRa(arr(2)), mRo(arr(3))
   {
   }

   Array EkQRaRoParameters::params() const
   {
      Array values(this->nParams());

      values(0) = this->E();
      values(1) = this->q();
      values(2) = this->Ra();
      values(3) = this->Ro();

      return values;
   }

   EPMFloat EkQRaRoParameters::nsDt() const
   {
      return this->Ro();
   }

   EPMFloat EkQRaRoParameters::nsDiffusion() const
   {
      return this->E();
   }

   EPMFloat EkQRaRoParameters::nsLorentz() const
   {
      return 1.0;
   }

   EPMFloat EkQRaRoParameters::nsAdvection() const
   {
      return this->Ro();
   }

   EPMFloat EkQRaRoParameters::nsBuoyancy() const
   {
      return this->q()*this->Ra();
   }

   EPMFloat EkQRaRoParameters::nsCoriolis() const
   {
      return 1.0;
   }

   EPMFloat EkQRaRoParameters::indDt() const
   {
      return 1.0;
   }

   EPMFloat EkQRaRoParameters::indDiffusion() const
   {
      return 1.0;
   }

   EPMFloat EkQRaRoParameters::indAdvection() const
   {
      return 1.0;
   }

   EPMFloat EkQRaRoParameters::tptDt() const
   {
      return 1.0;
   }

   EPMFloat EkQRaRoParameters::tptDiffusion() const
   {
      return this->q();
   }

   EPMFloat EkQRaRoParameters::tptAdvection() const
   {
      return 1.0;
   }

   EPMFloat EkQRaRoParameters::alfvenFactor() const
   {
      return this->Ro();
   }

   EPMFloat EkQRaRoParameters::alfvenDamping(const EPMFloat delta) const
   {
      return (this->E()+this->Ro())/(2.0*delta);
   }

   EPMFloat EkQRaRoParameters::meFactor() const
   {
      return 1.0/this->Ro();
   }

   void EkQRaRoParameters::testGlobalCFL(EPMFloat &rDt) const
   {
      if(this->Ro() != 0.0)
      {
         rDt = std::min(rDt, this->Ro());
      }

      if(this->E() != 0.0)
      {
         rDt = std::min(rDt, std::sqrt(this->E()));
      }
   }

}
