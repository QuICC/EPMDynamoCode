/** \file ETDOperators.hpp
 *  \brief General implemenation of a set of ETDOperators
 */

#ifndef ETDOPERATORS_HPP
#define ETDOPERATORS_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "General/EPMTypedefs.hpp"
#include "Operators/BoundedOperatorSet.hpp"

namespace EPMDynamo {

   /**
    * \brief General implemenation of a set of ETDOperators
    *
    * \tparam TSimType Type of the simulation
    * \tparam TOpType Type of the operator
    */
   template <typename TSimType, typename TOpType> class ETDOperators: public BoundedOperatorSet<TOpType>
   {
      public:
         /**
          * @brief Constructor
          *
          * @param pTrunc Truncation information
          */
         ETDOperators(SmartTruncation pTrunc);

         /**
          * @brief Destructor
          */
         virtual ~ETDOperators() {};
         
      protected:

      private:
   };

   template <typename TSimType> ETDOperators<TSimType>::ETDOperators(SmartTruncation pTrunc)
      : BoundedOperatorSet<TOpType>(pSTrunc)
   {
   }

}

#endif // ETDOPERATORS_HPP
