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
          * @param hasL0 Is l=0 mode required?
          */
         ETDOperators(SmartTruncation pTrunc, bool hasL0);

         /**
          * @brief Destructor
          */
         virtual ~ETDOperators() {};
         
      protected:

      private:
   };

   template <typename TSimType, typename TOpType> ETDOperators<TSimType, TOpType>::ETDOperators(SmartTruncation pTrunc, bool hasL0)
      : BoundedOperatorSet<TOpType>(pTrunc, hasL0)
   {
   }

}

#endif // ETDOPERATORS_HPP
