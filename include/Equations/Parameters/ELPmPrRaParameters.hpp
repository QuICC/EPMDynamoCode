/** \file ELPmPrRaParameters.hpp
 *  \brief Definition of the non-dimensionalisation based on the Ekman, Elsasser, magnetic Prandtl, Prandtl and Rayleigh numbers
 */

#ifndef ELPMPRRAPARAMETERS_HPP
#define ELPMPRRAPARAMETERS_HPP

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
    * \brief Definition of the non-dimensionalisation based on the Ekman, Elsasser, Prandtl, magnetic Prandtl and Rayleigh numbers
    */
   class ELPmPrRaParameters: public EquationParameters
   {
      public:
         /// Name (ID) of the equation parameters
         static const std::string  id;

         /**
         * @brief Constructor
         *
         * \param E Ekman number
         * \param L Elsasser number
         * \param Pm magnetic Prandtl number
         * \param Pr Prandtl number
         * \param Ra Rayleigh number
         */
         ELPmPrRaParameters(EPMFloat E, EPMFloat L, EPMFloat Pm, EPMFloat Pr, EPMFloat Ra);

         /**
         * @brief Constructor
         *
         * \param arr Array of non-dimensional parameters: (0) Ekman, (1) Elsasser, (2) magnetic Prandtl, (3) Prandtl, (4) Rayleigh
         */
         ELPmPrRaParameters(const Array& arr);

         /**
         * @brief Destructor
         */
         virtual ~ELPmPrRaParameters() {};

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
          * @brief Get the Elsasser number
          */
         EPMFloat L() const;

         /**
          * @brief Get the magnetic Prandtl number
          */
         EPMFloat Pm() const; 

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
          * @brief Get the magnetic energy normalisation factor
          */
         EPMFloat meFactor() const;

         /**
          * @brief Test the given timestep value against global CFL conditions
          */
         void testGlobalCFL(EPMFloat &rDt) const;

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
          * @brief The Elsasser number
          */
         EPMFloat mL;

         /**
          * @brief The magnetic Prandtl number
          */
         EPMFloat mPm;

         /**
          * @brief The Prandtl number
          */
         EPMFloat mPr;

         /**
          * @brief The Rayleigh number
          */
         EPMFloat mRa;
   };

   inline std::string ELPmPrRaParameters::type() const
   {
      return ELPmPrRaParameters::id;
   }

   inline int ELPmPrRaParameters::nParams() const
   {
      return 5;
   }

   inline EPMFloat ELPmPrRaParameters::E() const
   {
      return this->mE;
   }

   inline EPMFloat ELPmPrRaParameters::L() const
   {
      return this->mL;
   }

   inline EPMFloat ELPmPrRaParameters::Pm() const
   {
      return this->mPm;
   }

   inline EPMFloat ELPmPrRaParameters::Pr() const
   {
      return this->mPr;
   }

   inline EPMFloat ELPmPrRaParameters::Ra() const
   {
      return this->mRa;
   }
}

#endif // ELPMPRRAPARAMETERS_HPP
