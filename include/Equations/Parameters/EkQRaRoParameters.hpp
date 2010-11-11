/** \file EkQRaRoParameters.hpp
 *  \brief Definition of the non-dimensionalisation based on the Ekman, Roberts, Rayleigh and magnetic Rossby numbers
 */

#ifndef EKQRAROPARAMETERS_HPP
#define EKQRAROPARAMETERS_HPP

// Configuration includes
//
#include "Config/SmartPointer.h"

// System includes
//

// External includes
//

// Project includes
//
#include "General/EPMTypedefs.hpp"
#include "Equations/Parameters/EquationParameters.hpp"

namespace EPMDynamo {

   /**
    * @brief This class provides the constant coefficients apparearing in the equations
    *
    * It could also later on, provide conversion methods between non-dimensional values and 
    * real physical values.
    */
   class EkQRaRoParameters: public EquationParameters
   {
      public:
         /**
         * @brief Constructor
         *
         * \param E Ekman number
         * \param q Roberts number
         * \param Ra Rayleigh number
         * \param Ro Rossby number
         */
         EkQRaRoParameters(EPMFloat E, EPMFloat q, EPMFloat Ra, EPMFloat Ro);

         /**
         * @brief Constructor
         *
         * \param arr Array of non-dimensional parameters: (0) Ekman, (1) Roberts, (2) Rayleigh, (3) Rossby
         */
         EkQRaRoParameters(const Array& arr);

         /**
         * @brief Destructor
         */
         virtual ~EkQRaRoParameters() {};

         /**
          * @brief Type of the parameters
          */
         std::string type() const;

         /**
          * @brief Get the number of parameters
          */
         int nParams() const;

         /**
          * @brief Get the values of all the parameters
          */
         Array params() const;

         /**
          * @brief Get the Rossby number
          */
         EPMFloat Ro() const; 

         /**
          * @brief Get the Ekman number
          */
         EPMFloat E() const;

         /**
          * @brief Get the q number
          */
         EPMFloat q() const;

         /**
          * @brief Get the Rayleigh number
          */
         EPMFloat Ra() const;

         /**
          * @brief Get non dimensionalisation parameter for time derivative in Navier-Stokes equation
          */
         EPMFloat nsDt() const;

         /**
          * @brief Get non dimensionalisation parameter for diffusion term in Navier-Stokes equation
          */
         EPMFloat nsDiffusion() const;

         /**
          * @brief Get non dimensionalisation parameter for Lorentz force term in Navier-Stokes equation
          */
         EPMFloat nsLorentz() const;

         /**
          * @brief Get non dimensionalisation parameter for advection term in Navier-Stokes equation
          */
         EPMFloat nsAdvection() const;

         /**
          * @brief Get non dimensionalisation parameter for buoyancy term in Navier-Stokes equation
          */
         EPMFloat nsBuoyancy() const;

         /**
          * @brief Get non dimensionalisation parameter for Coriolis term in Navier-Stokes equation
          */
         EPMFloat nsCoriolis() const;

         /**
          * @brief Get non dimensionalisation parameter for time derivative in induction equation
          */
         EPMFloat indDt() const;

         /**
          * @brief Get non dimensionalisation parameter for diffusion tern in induction equation
          */
         EPMFloat indDiffusion() const;

         /**
          * @brief Get non dimensionalisation parameter for advection term in induction equation
          */
         EPMFloat indAdvection() const;

         /**
          * @brief Get non dimensionalisation parameter for time derivative in transport equation
          */
         EPMFloat tptDt() const;

         /**
          * @brief Get non dimensionalisation parameter for diffusion tern in transport equation
          */
         EPMFloat tptDiffusion() const;

         /**
          * @brief Get non dimensionalisation parameter for advection term in transport equation
          */
         EPMFloat tptAdvection() const;

         /**
          * @brief Damping factor of the Alfven velocity used in CFL condition
          * \f[
          *    U_R = B^2/\sqrt{factor*B^2 + damping}
          * \f]
          *
          * @param delta Smallest physical scale
          */
         EPMFloat alfvenDamping(const EPMFloat delta) const;

         /**
          * @brief Alfven velocity multiplicative factor. The velocity is written in following form 
          * \f[
          *    U_R = B^2/\sqrt{factor*B^2 + damping}
          * \f]
          */
         EPMFloat alfvenFactor() const;

         /**
          * @brief Get the magnetic energy normalisation factor
          */
         EPMFloat meFactor() const;

         /**
          * @brief Test the given timestep value against global CFL conditions
          */
         void testGlobalCFL(EPMFloat &rDt) const;
         
      protected:

      private:

         /**
          * @brief The Ekman number
          */
         EPMFloat mE;

         /**
          * @brief The Roberts q number
          */
         EPMFloat mQ;

         /**
          * @brief The Rayleigh number
          */
         EPMFloat mRa;

         /**
          * @brief The Rossby number
          */
         EPMFloat mRo;
   };

   inline std::string EkQRaRoParameters::type() const
   {
      return "EkQRaRo";
   }

   inline int EkQRaRoParameters::nParams() const
   {
      return 4;
   }

   inline EPMFloat EkQRaRoParameters::Ro() const
   {
      return this->mRo;
   }

   inline EPMFloat EkQRaRoParameters::E() const
   {
      return this->mE;
   }

   inline EPMFloat EkQRaRoParameters::q() const
   {
      return this->mQ;
   }

   inline EPMFloat EkQRaRoParameters::Ra() const
   {
      return this->mRa;
   }
}

#endif // EKQRAROPARAMETERS_HPP
