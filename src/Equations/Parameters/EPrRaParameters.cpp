/** \file EPrRaParameters.cpp
 *  \brief Source of the implementation of the equation non dimensional parameters
 */

// System includes
//

// External includes
//

// Class include
//
#include "Equations/Parameters/EPrRaParameters.hpp"

// Project includes
//

namespace EPMDynamo {

   const std::string EPrRaParameters::id = "EPrRa";

   EPrRaParameters::EPrRaParameters(EPMFloat E, EPMFloat Pr, EPMFloat Ra)
      : EquationParameters(), mE(E), mPr(Pr), mRa(Ra)
   {
   }

   EPrRaParameters::EPrRaParameters(const Array& arr)
      : EquationParameters(), mE(arr(0)), mPr(arr(1)), mRa(arr(2))
   {
   }

   Array EPrRaParameters::params() const
   {
      Array values(this->nParams());

      values(0) = this->E();
      values(1) = this->Pr();
      values(2) = this->Ra();

      return values;
   }

   EPMFloat EPrRaParameters::nsDt() const
   {
      return this->E();
   }

   EPMFloat EPrRaParameters::nsDiffusion() const
   {
      return this->E();
   }

   EPMFloat EPrRaParameters::nsLorentz() const
   {
      return 0.0;
   }

   EPMFloat EPrRaParameters::nsAdvection() const
   {
      //return 1.0;
      //Leo: test
      return this->E();
   }

   EPMFloat EPrRaParameters::nsBuoyancy() const
   {  
      //Leo: standard Rayleigh
      //return this->E()*this->Ra();
      //Leo: modified Rayleigh
      return this->Ra();
   }

   EPMFloat EPrRaParameters::nsCoriolis() const
   {
      return 1.0;
   }

   EPMFloat EPrRaParameters::indDt() const
   {
      return 0.0;
   }

   EPMFloat EPrRaParameters::indDiffusion() const
   {
      return 0.0;
   }

   EPMFloat EPrRaParameters::indAdvection() const
   {
      return 0.0;
   }

   EPMFloat EPrRaParameters::tptDt() const
   {
      return this->Pr();
   }

   EPMFloat EPrRaParameters::tptDiffusion() const
   {
      return 1.0;
   }

   EPMFloat EPrRaParameters::tptAdvection() const
   {
      //return 1.0;
	  return this->Pr();
   }

   EPMFloat EPrRaParameters::alfvenFactor() const
   {
      return 0.0;
   }

   EPMFloat EPrRaParameters::alfvenDamping(const EPMFloat delta) const
   {
      return 0.0;
   }

   EPMFloat EPrRaParameters::keFactor() const
   {
      return 1.0;
   }

   EPMFloat EPrRaParameters::meFactor() const
   {
      return 0.0;
   }

   EPMFloat EPrRaParameters::inertialCFL() const
   {
      //return 1.0;
      //Leo
      //Inertial waves
      //EParameters.cpp
      //1/(2*Omega)
      //Delta Tg= C/10*1/(2*Omega)
      //Precession: Time=1/Omega
      //return 0.5/10.;
      //Convection: Time = r^2/nu
      return this->E()/10.;
   }

   EPMFloat EPrRaParameters::torsionalCFL() const
   {
      //For magnetic
      return 1.0;
   }

   EPMFloat EPrRaParameters::codSourceScale() const
   {
      return 1.0;
   }

   EPMFloat EPrRaParameters::imposedMagScale() const
   {
      return 0.0;
   }

}
