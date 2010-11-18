/** \file EquationParameters.hpp
 *  \brief Definition of the non-dimensional parameters
 */

#ifndef EQUATIONPARAMETERS_HPP
#define EQUATIONPARAMETERS_HPP

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

namespace EPMDynamo {

   /**
    * @brief This class provides the constant coefficients apparearing in the equations
    */
   class EquationParameters
   {
      public:
         /**
         * @brief Constructor
         */
         EquationParameters();

         /**
         * @brief Destructor
         */
         virtual ~EquationParameters() {};

         /**
          * @brief Type of the parameters
          */
         virtual std::string type() const = 0;

         /**
          * @brief Get the number of parameters
          */
         virtual int nParams() const = 0;

         /**
          * @brief Get the values of all the parameters
          */
         virtual Array params() const = 0;

         /**
          * @brief Get non dimensionalisation parameter for time derivative in Navier-Stokes equation
          */
         virtual EPMFloat nsDt() const = 0;

         /**
          * @brief Get non dimensionalisation parameter for diffusion term in Navier-Stokes equation
          */
         virtual EPMFloat nsDiffusion() const = 0;

         /**
          * @brief Get non dimensionalisation parameter for Lorentz force term in Navier-Stokes equation
          */
         virtual EPMFloat nsLorentz() const = 0;

         /**
          * @brief Get non dimensionalisation parameter for advection term in Navier-Stokes equation
          */
         virtual EPMFloat nsAdvection() const = 0;

         /**
          * @brief Get non dimensionalisation parameter for buoyancy term in Navier-Stokes equation
          */
         virtual EPMFloat nsBuoyancy() const = 0;

         /**
          * @brief Get non dimensionalisation parameter for Coriolis term in Navier-Stokes equation
          */
         virtual EPMFloat nsCoriolis() const = 0;

         /**
          * @brief Get non dimensionalisation parameter for time derivative in induction equation
          */
         virtual EPMFloat indDt() const = 0;

         /**
          * @brief Get non dimensionalisation parameter for diffusion tern in induction equation
          */
         virtual EPMFloat indDiffusion() const = 0;

         /**
          * @brief Get non dimensionalisation parameter for advection term in induction equation
          */
         virtual EPMFloat indAdvection() const = 0;

         /**
          * @brief Get non dimensionalisation parameter for time derivative in transport equation
          */
         virtual EPMFloat tptDt() const = 0;

         /**
          * @brief Get non dimensionalisation parameter for diffusion tern in transport equation
          */
         virtual EPMFloat tptDiffusion() const = 0;

         /**
          * @brief Get non dimensionalisation parameter for advection term in transport equation
          */
         virtual EPMFloat tptAdvection() const = 0;

         /**
          * @brief Alfven velocity multiplicative factor. The velocity is written in following form 
          * \f[
          *    U_R = B^2/\sqrt{factor*B^2 + damping}
          * \f]
          */
         virtual EPMFloat alfvenFactor() const = 0;

         /**
          * @brief Damping factor of the Alfven velocity used in CFL condition
          * \f[
          *    U_R = B^2/\sqrt{factor*B^2 + damping}
          * \f]
          *
          * @param delta Smallest physical scale
          */
         virtual EPMFloat alfvenDamping(const EPMFloat delta) const = 0;

         /**
          * @brief Get the magnetic energy normalisation factor
          */
         virtual EPMFloat meFactor() const = 0;

         /**
          * @brief Test the given timestep value against global CFL conditions
          */
         virtual void testGlobalCFL(EPMFloat &rDt) const = 0;

         /**
          * @brief Scaling factor for the codensity source field
          */
         virtual EPMFloat codSourceScale() const = 0;

         /**
          * @brief Scaling factor for the imposed magnetic field
          */
         virtual EPMFloat imposedMagScale() const = 0;
         
      protected:

      private:
   };

   /// Typedef for geting a singleton of the EquationParameters
   typedef EPMSHARED_PTR<EquationParameters>   SmartEqParameters;
}

#endif // EQUATIONPARAMETERS_HPP
