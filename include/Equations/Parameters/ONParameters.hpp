/** \file ONParameters.hpp
 *  \brief Definition of the non-dimensionalisation based on Omega and nu
 */

#ifndef ONPARAMETERS_HPP
#define ONPARAMETERS_HPP

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
    * \brief Definition of the non-dimensionalisation based on the Omega and nu
    */
   class ONParameters: public EquationParameters
   {
      public:
         /// Name (ID) of the equation parameters
         static const std::string  id;

         /**
         * @brief Constructor
         *
         * \param Omega omega
         * \param nu viscosity
         */
         ONParameters(EPMFloat omega, EPMFloat nu);

         /**
         * @brief Constructor
         *
         * \param arr Array of non-dimensional parameters: (0) Ekman, (1) magnetic Prandtl
         */
         ONParameters(const Array& arr);

         /**
         * @brief Destructor
         */
         virtual ~ONParameters() {};

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
          * @brief Get the Omega
          */
         EPMFloat omega() const;

         /**
          * @brief Get the viscosity
          */
         EPMFloat nu() const; 

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
          * @brief The Omega
          */
         EPMFloat mOmega;

         /**
          * @brief The viscosity
          */
         EPMFloat mNu;
   };

   inline std::string ONParameters::type() const
   {
      return ONParameters::id;
   }

   inline int ONParameters::nParams() const
   {
      return 2;
   }

   inline EPMFloat ONParameters::omega() const
   {
      return this->mOmega;
   }

   inline EPMFloat ONParameters::nu() const
   {
      return this->mNu;
   }
}

#endif // ONPARAMETERS_HPP
