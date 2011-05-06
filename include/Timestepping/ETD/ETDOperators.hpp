/** \file ETDOperators.hpp
 *  \brief General implementation of a set of ETDOperators
 */

#ifndef ETDOPERATORS_HPP
#define ETDOPERATORS_HPP

// Configuration includes
//

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
    * \tparam TOpType Type of the operator
    */
   template <typename TOpType> class ETDOperators: public BoundedOperatorSet<TOpType>
   {
      public:
         /**
          * @brief Constructor
          *
          * @param c Stiffness constant
          * @param pTrunc Truncation information
          * @param hasL0 Is l=0 mode required?
          */
         ETDOperators(EPMFloat c, SmartTruncation pTrunc, bool hasL0);

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
          * @brief Storage for the stiffness constant c
          */
         EPMFloat mC;

         /**
          * @brief Storage for the current timestep length
          */
         EPMFloat mH;
   };

   template <typename TOpType> inline EPMFloat ETDOperators<TOpType>::h() const
   {
      // Guard from uninitialised timestep length
      assert(this->mH > 0.0);

      return this->mH;
   }

   template <typename TOpType> ETDOperators<TOpType>::ETDOperators(EPMFloat c, SmartTruncation pTrunc, bool hasL0)
      : BoundedOperatorSet<TOpType>(pTrunc, hasL0), mC(c), mH(-1)
   {
   }

   template <typename TOpType> void ETDOperators<TOpType>::updateTimestep(EPMFloat h)
   {
      // Guard from uninitialised timestep length
      assert(h > 0.0);

      this->mH = h;
   }

}

#endif // ETDOPERATORS_HPP
