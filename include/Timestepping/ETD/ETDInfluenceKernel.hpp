/** \file ETDInfluenceKernel.hpp
 *  \brief Implemenation of influence of the kernel solution
 */

#ifndef ETDINFLUENCEKERNEL_HPP
#define ETDINFLUENCEKERNEL_HPP

// Configuration includes
//
#include "Config/SimulationConfig.hpp"

// System includes
//

// External includes
//

// Project includes
//
#include "Timestepping/ETD/ETDSchemeTraits.hpp"

namespace EPMDynamo {

   /**
    * \brief Implemenation of influence of the kernel solution
    */
   class ETDInfluenceKernel
   {
      public:
         /// Typedef from Simulation trait to local scalar type
         typedef SimulationConfig::NumericalScheme::ScalarType    ScalarType;

         /// Typedef for a smart pointer to ETDOperators
         typedef EPMSHARED_PTR<ETDSchemeTraits::Operators> SmartETDOperators;

         /**
          * @brief Constructor
          *
          * @param nFactor Non linear terms multiplicative factor
          * @param pOpM1 Pointer to operator M1
          * @param pOpM2 Pointer to operator M2
          */
         ETDInfluenceKernel(EPMFloat nFactor, SmartETDOperators pOpM1);

         /**
          * @brief Destructor
          */
         virtual ~ETDInfluenceKernel() {};

         /**
          * @brief Compute the influence of the kernel
          *
          * @param rKernel Kernel decomposition
          * @param l Harmonic degree
          */
         void computeInfluence(Array  &rKernel, const int l);
         
      protected:
         /**
          * @brief multiplicative factor required for the non linear terms
          */
         EPMFloat mNFactor;

         /**
          * @brief Pointer to the ETD operator M1
          */
         SmartETDOperators  mpOpM1;

      private:
   };

}

#endif // ETDINFLUENCEKERNEL_HPP
