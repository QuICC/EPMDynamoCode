/** \file EPmQRaParameters.hpp
 *  \brief Definition of the non-dimensionalisation based on the Ekman, magnetic Prandtl, Roberts and Rayleigh numbers
 */

#ifndef EPMQRAPARAMETERS_HPP
#define EPMQRAPARAMETERS_HPP

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
    * \brief Definition of the non-dimensionalisation based on the Ekman, magnetic Prandtl, Roberts and Rayleigh numbers
    */
   class EPmQRaParameters: public EquationParameters
   {
      public:
         /// Name (ID) of the equation parameters
         static const std::string  id;

         /**
         * @brief Constructor
         *
         * \param E Ekman number
         * \param Pm magnetic Prandtl number
         * \param q Roberts number
         * \param Ra Rayleigh number
         */
         EPmQRaParameters(EPMFloat E, EPMFloat Pm, EPMFloat q, EPMFloat Ra);

         /**
         * @brief Constructor
         *
         * \param arr Array of non-dimensional parameters: (0) Ekman, (1) magnetic Prandtl, (2) Roberts, (3) Rayleigh
         */
         EPmQRaParameters(const Array& arr);

         /**
         * @brief Destructor
         */
         virtual ~EPmQRaParameters() {};

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
          * @brief Get the Ekman number
          */
         EPMFloat E() const;

         /**
          * @brief Get the magnetic Prandtl number
          */
         EPMFloat Pm() const; 

         /**
          * @brief Get the Roberts number
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
          * @brief Get the kinetic energy normalisation factor
          */
         EPMFloat keFactor() const;

         /**
          * @brief Get the magnetic energy normalisation factor
          */
         EPMFloat meFactor() const;

         /**
          * @brief Get the global inertial wave CFL condition
          */
         EPMFloat inertialCFL() const;

         /**
          * @brief Get the global torsional osciallations CFL condition
          */
         EPMFloat torsionalCFL() const;

         /**
          * @brief Scaling factor for the codensity source field
          */
         EPMFloat codSourceScale() const;

         /**
          * @brief Scaling factor for the imposed magnetic field
          */
         EPMFloat imposedMagScale() const;
         
      protected:

      private:

         /**
          * @brief The Ekman number
          */
         EPMFloat mE;

         /**
          * @brief The magnetic Prandtl number
          */
         EPMFloat mPm;

         /**
          * @brief The Roberts number
          */
         EPMFloat mQ;

         /**
          * @brief The Rayleigh number
          */
         EPMFloat mRa;
   };

   inline std::string EPmQRaParameters::type() const
   {
      return EPmQRaParameters::id;
   }

   inline int EPmQRaParameters::nParams() const
   {
      return 4;
   }

   inline EPMFloat EPmQRaParameters::E() const
   {
      return this->mE;
   }

   inline EPMFloat EPmQRaParameters::Pm() const
   {
      return this->mPm;
   }

   inline EPMFloat EPmQRaParameters::q() const
   {
      return this->mQ;
   }

   inline EPMFloat EPmQRaParameters::Ra() const
   {
      return this->mRa;
   }
}

#endif // EPMQRAPARAMETERS_HPP
