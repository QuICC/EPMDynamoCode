/** \file EPrRaParameters.hpp
    * \brief Definition of the non-dimensionalisation for rotating thermal convection. Based on Jones 2000 paper. It needs the Ekman number. the Prandtl number Pr and a rayleigh number
 */

#ifndef EPRRAPARAMETERS_HPP
#define EPRRAPARAMETERS_HPP

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
    * \brief Definition of the non-dimensionalisation for rotating thermal convection. Based on Jones 2000 paper. It needs the Ekman number. the Prandtl number Pr and a rayleigh number
    */
   class EPrRaParameters: public EquationParameters
   {
      public:
         /// Name (ID) of the equation parameters
         static const std::string  id;

         /**
         * @brief Constructor
         *
         * \param E Ekman number
         * \param Pr Prandtl number
         * \param Ra Rayleigh number
         */
         EPrRaParameters(EPMFloat E, EPMFloat Pr, EPMFloat Ra);

         /**
         * @brief Constructor
         *
         * \param arr Array of non-dimensional parameters: (0) Ekman, (1) Prandtl, (2) Rayleigh
         */
         EPrRaParameters(const Array& arr);

         /**
         * @brief Destructor
         */
         virtual ~EPrRaParameters() {};

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
          * @brief Get the Prandtl number
          */
         EPMFloat Pr() const; 

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
          * @brief The Prandtl number
          */
         EPMFloat mPr;

         /**
          * @brief The Rayleigh number
          */
         EPMFloat mRa;
   };

   inline std::string EPrRaParameters::type() const
   {
      return EPrRaParameters::id;
   }

   inline int EPrRaParameters::nParams() const
   {
      return 3;
   }

   inline EPMFloat EPrRaParameters::E() const
   {
      return this->mE;
   }

   inline EPMFloat EPrRaParameters::Pr() const
   {
      return this->mPr;
   }

   inline EPMFloat EPrRaParameters::Ra() const
   {
      return this->mRa;
   }
}

#endif // EPRRAPARAMETERS_HPP
