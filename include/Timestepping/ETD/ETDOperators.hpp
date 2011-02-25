/** \file ETDOperators.hpp
 *  \brief General implementation of a set of ETDOperators
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
    * \brief General implementation of a set of ETDOperators
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

         /**
          * @brief Get the current timestep length
          */
         EPMFloat h() const;

         /**
          * @brief Update the timestep length
          *
          * @param h New timestep length
          */
         void updateTimestep(const EPMFloat h);
         
      protected:

      private:
         /**
          * @brief Storage for the current timestep length
          */
         EPMFloat mH;
   };

   template <typename TSimType, typename TOpType> inline EPMFloat ETDOperators<TSimType, TOpType>::h() const
   {
      // Guard from uninitialised timestep length
      assert(this->mH > 0.0);

      return this->mH;
   }

   template <typename TSimType, typename TOpType> ETDOperators<TSimType, TOpType>::ETDOperators(SmartTruncation pTrunc, bool hasL0)
      : BoundedOperatorSet<TOpType>(pTrunc, hasL0), mH(-1)
   {
   }

   template <typename TSimType, typename TOpType> void ETDOperators<TSimType, TOpType>::updateTimestep(EPMFloat h)
   {
      this->mH = h;
   }

}

#endif // ETDOPERATORS_HPP
