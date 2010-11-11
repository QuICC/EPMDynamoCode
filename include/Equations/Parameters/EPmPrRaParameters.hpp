/** \file EPmPrRaParameters.hpp
 *  \brief Definition of the non-dimensionalisation based on the Ekman, magnetic Prandtl, Prandtl and Rayleigh
 */

#ifndef EPMPRRAPARAMETERS_HPP
#define EPMPRRAPARAMETERS_HPP

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
    * \brief Definition of the non-dimensionalisation based on the Ekman, Prandtl, magnetic Prandtl and Rayleigh
    */
   class EPmPrRaParameters: public EquationParameters
   {
      public:
         /// Name (ID) of the equation parameters
         static const std::string  id;

         /**
         * @brief Constructor
         *
         * \param E Ekman number
         * \param Pm magnetic Prandtl number
         * \param Pr Prandtl number
         * \param Ra Rayleigh number
         */
         EPmPrRaParameters(EPMFloat E, EPMFloat Pm, EPMFloat Pr, EPMFloat Ra);

         /**
         * @brief Constructor
         *
         * \param arr Array of non-dimensional parameters: (0) Ekman, (1) magnetic Prandtl, (2) Prandtl, (3) Rayleigh
         */
         EPmPrRaParameters(const Array& arr);

         /**
         * @brief Destructor
         */
         virtual ~EPmPrRaParameters() {};

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
          * @brief The Prandtl number
          */
         EPMFloat mPr;

         /**
          * @brief The Rayleigh number
          */
         EPMFloat mRa;
   };

   inline std::string EPmPrRaParameters::type() const
   {
      return EPmPrRaParameters::id;
   }

   inline int EPmPrRaParameters::nParams() const
   {
      return 4;
   }

   inline EPMFloat EPmPrRaParameters::E() const
   {
      return this->mE;
   }

   inline EPMFloat EPmPrRaParameters::Pm() const
   {
      return this->mPm;
   }

   inline EPMFloat EPmPrRaParameters::Pr() const
   {
      return this->mPr;
   }

   inline EPMFloat EPmPrRaParameters::Ra() const
   {
      return this->mRa;
   }
}

#endif // EPMPRRAPARAMETERS_HPP
