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
    *
    * It could also later on, provide conversion methods between non-dimensional values and 
    * real physical values.
    */
   class EquationParameters
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
         EquationParameters(EPMFloat E, EPMFloat q, EPMFloat Ra, EPMFloat Ro);

         /**
         * @brief Constructor
         *
         * \param arr Array of non-dimensional parameters: (0) Ekman, (1) Roberts, (2) Rayleigh, (3) Rossby
         */
         EquationParameters(const Array& arr);

         /**
         * @brief Destructor
         */
         virtual ~EquationParameters() {};

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

   inline EPMFloat EquationParameters::Ro() const
   {
      return this->mRo;
   }

   inline EPMFloat EquationParameters::E() const
   {
      return this->mE;
   }

   inline EPMFloat EquationParameters::q() const
   {
      return this->mQ;
   }

   inline EPMFloat EquationParameters::Ra() const
   {
      return this->mRa;
   }

   /// Typedef for geting a singleton of the EquationParameters
   typedef EPMSHARED_PTR<EquationParameters>   SmartEqParameters;
}

#endif // EQUATIONPARAMETERS_HPP
