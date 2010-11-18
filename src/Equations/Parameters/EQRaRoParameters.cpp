/** \file EQRaRoParameters.cpp
 *  \brief Source of the implementation of the equation non dimensional parameters
 */

// System includes
//

// External includes
//

// Class include
//
#include "Equations/Parameters/EQRaRoParameters.hpp"

// Project includes
//

namespace EPMDynamo {

   const std::string EQRaRoParameters::id = "EQRaRo";

   EQRaRoParameters::EQRaRoParameters(EPMFloat E, EPMFloat q, EPMFloat Ra, EPMFloat Ro)
      : EquationParameters(), mE(E), mQ(q), mRa(Ra), mRo(Ro)
   {
   }

   EQRaRoParameters::EQRaRoParameters(const Array& arr)
      : EquationParameters(), mE(arr(0)), mQ(arr(1)), mRa(arr(2)), mRo(arr(3))
   {
   }

   Array EQRaRoParameters::params() const
   {
      Array values(this->nParams());

      values(0) = this->E();
      values(1) = this->q();
      values(2) = this->Ra();
      values(3) = this->Ro();

      return values;
   }

   EPMFloat EQRaRoParameters::nsDt() const
   {
      return this->Ro();
   }

   EPMFloat EQRaRoParameters::nsDiffusion() const
   {
      return this->E();
   }

   EPMFloat EQRaRoParameters::nsLorentz() const
   {
      return 1.0;
   }

   EPMFloat EQRaRoParameters::nsAdvection() const
   {
      return this->Ro();
   }

   EPMFloat EQRaRoParameters::nsBuoyancy() const
   {
      return this->q()*this->Ra();
   }

   EPMFloat EQRaRoParameters::nsCoriolis() const
   {
      return 1.0;
   }

   EPMFloat EQRaRoParameters::indDt() const
   {
      return 1.0;
   }

   EPMFloat EQRaRoParameters::indDiffusion() const
   {
      return 1.0;
   }

   EPMFloat EQRaRoParameters::indAdvection() const
   {
      return 1.0;
   }

   EPMFloat EQRaRoParameters::tptDt() const
   {
      return 1.0;
   }

   EPMFloat EQRaRoParameters::tptDiffusion() const
   {
      return this->q();
   }

   EPMFloat EQRaRoParameters::tptAdvection() const
   {
      return 1.0;
   }

   EPMFloat EQRaRoParameters::alfvenFactor() const
   {
      return this->Ro();
   }

   EPMFloat EQRaRoParameters::alfvenDamping(const EPMFloat delta) const
   {
      return (this->E()+this->Ro())/(2.0*delta);
   }

   EPMFloat EQRaRoParameters::meFactor() const
   {
      return 1.0/this->Ro();
   }

   void EQRaRoParameters::testGlobalCFL(EPMFloat &rDt) const
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

   EPMFloat EQRaRoParameters::codSourceScale() const
   {
      return this->q();
   }

   EPMFloat EQRaRoParameters::imposedMagScale() const
   {
      return 1.0;
   }

}
