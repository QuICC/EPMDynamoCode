/** \file Integrator.hpp
 *  \brief Implementation of a initial state generator
 */

#ifndef INTEGRATOR_HPP
#define INTEGRATOR_HPP

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
#include "Utilities/GeneratorBase.hpp"

namespace EPMDynamo {

   /**
    * @brief This class implements a few methods for easier implementation of an initial state generator
    *
    * \tparam TSimTraits Traits for the generator
    */
   template <typename TSimTraits> class Integrator: public GeneratorBase<TSimTraits>
   {
      public:
         /**
          * @brief Constructor
          */
         Integrator();

         /**
          * @brief Simple empty destructor
          */
         virtual ~Integrator() {};

         /**
          * @brief Integrate field
          */
         void integrate();

      protected:


      private:
   };

   template <typename TSimTraits> Integrator<TSimTraits>::Integrator()
      : GeneratorBase<TSimTraits>()
   {
   }

   template <typename TSimTraits> void Integrator<TSimTraits>::integrate()
   {
      if(TSimTraits::NeedCodensity)
      {
         this->codC().rOc().updateSpectra();
         std::cerr << "Codensity: " << this->codC().oc().energy().transpose() << std::endl;
      }

      if(TSimTraits::NeedVelocity)
      {
         this->velV().rOc().updateSpectra();
         std::cerr << "Velocity: " << this->velV().oc().energy().transpose() << std::endl;
      }

      if(TSimTraits::NeedMagnetic)
      {
         this->magB().rOc().updateSpectra();
         std::cerr << "Magnetic: " << this->magB().oc().energy().transpose() << std::endl;
      }
   }

}

#endif // INTEGRATOR_HPP
