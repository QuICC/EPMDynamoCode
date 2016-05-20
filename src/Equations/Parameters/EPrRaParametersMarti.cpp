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
	  //Leo:
	  //Thermal Convection Benchmark 1
	  //Thermal Precession
      //Time: nu/(omega*d^2)
      return this->E();
      //Time: 1/omega
      //return 1.0;
   }

   EPMFloat EPrRaParameters::nsDiffusion() const
   {
	  //Leo:
	  //Thermal Convection Benchmark 1
	  //Thermal Precession
      //time: diffusion and rotation
      return this->E();
   }

   EPMFloat EPrRaParameters::nsLorentz() const
   {
      return 0.0;
   }

   EPMFloat EPrRaParameters::nsAdvection() const
   {
      //return 1.0;
      //Leo:
	  //Thermal Convection Benchmark 1
	  //Thermal Precession
      //Time: diffusion
      return this->E();
      //Time: rotation
      //return 1.0;
   }

   EPMFloat EPrRaParameters::nsBuoyancy() const
   {  
      //Leo: standard Rayleigh
      //return this->E()*this->Ra();
      //Leo:
	  //Thermal Convection Benchmark 1
	  //Thermal Precession
      //Time: diffusion
      //Modified Rayleigh (Marti)
      return this->Ra();
	  //Modified Rayleigh (Wei & Tilgner)
      //return this->Ra()/this->E();
      //Time: rotation 
      //return this->Ra();
   }

   EPMFloat EPrRaParameters::nsCoriolis() const
   {
	  //Leo: Thermal Convection Benchmark 1
      return 1.0;


	  //Thermal Precession
      //time:diffusion and rotation
	  //return 2.0;
 
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
	  //Leo:
	  //Thermal Convection Benchmark 1
	  //Thermal Precession
      return this->Pr();
   }

   EPMFloat EPrRaParameters::tptDiffusion() const
   {
	  //Leo
	  //Thermal Convection Benchmark 1
	  //Thermal Precession
      //time: diffusion
      return 1.0;
      //time: rotation
      //return this->E();
   }

   EPMFloat EPrRaParameters::tptAdvection() const
   {
      //return 1.0;
	  //Leo:
	  //Thermal Convection Benchmark 1
	  //Thermal Precession
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
   //Kinetic energy factor
   EPMFloat EPrRaParameters::keFactor() const
   {
	   return 0.5;
	   //return 1.0;
   }
   //Magnetic energy factor
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

	  //Thermal convection
      //Thermal-Convection: Time = r^2/nu
	  // T = 1/(2*Omega) = 1 / (2*Omega*r^2/nu) = E
      return this->E()/10.;

	  //Thermal Precession
      //time diffusion
	  //return this->E()/20.;
      //time rotation
      //return 0.5/10.;

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
