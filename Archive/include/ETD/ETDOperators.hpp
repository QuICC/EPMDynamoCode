/** \file ETDOperators.hpp
 */

#ifndef ETDOPERATORS_HPP
#define ETDOPERATORS_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "Simulations/SimulationTraits.hpp"
#include "General/Typedefs.hpp"
#include "Operators/HomogeneousSet.hpp"

namespace EPMDynamo {

   template <typename TSim> class ETDOperators: public HomogeneousSet<typename SimulationTraits<TSim>::OperatorType>
   {
      public:
         /**
          * @brief Constructor
          */
         ETDOperators(SmartSTrunc pTrunc);

         /**
          * @brief Destructor
          */
         virtual ~ETDOperators();
         
      protected:

      private:
   };

   template <typename TSim> ETDOperators<TSim>::ETDOperators(SmartSTrunc pTrunc)
      : HomogeneousSet<typename SimulationTraits<TSim>::OperatorType>(pSTrunc)
   {
   }

   template <typename TSim> ETDOperators<TSim>::~ETDOperators()
   {
   }

}

#endif // ETDOPERATORS_HPP
