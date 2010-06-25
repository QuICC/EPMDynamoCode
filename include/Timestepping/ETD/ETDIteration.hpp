/** \file ETDIteration.hpp
 *  \brief General implementation of an ETD iteration (step)
 */

#ifndef ETDITERATION_HPP
#define ETDITERATION_HPP

// System includes
//

// External includes
//

// Project includes
//

namespace EPMDynamo {

   /**
    * \brief General implementation of an ETD iteration (step)
    *
    * \tparam TSimType Type of the simulation
    */
   template <typename TSimType> class ETDIteration
   {
      public:
         /// Typedef from Simulation trait to local scalar type
         typedef typename TSimType::ScalarType    ScalarType;

         /**
          * @brief Constructor
          */
         ETDIteration();

         /**
          * @brief Destructor
          */
         virtual ~ETDIteration() {};

         /**
          * @brief Compute the implemented iteration
          *
          * @param rVar Input/Output variable
          * @param newNTerms New non linear terms
          */
         virtual void compute(ScalarType &rVar, ScalarType rNTerms) = 0;
         
      protected:

      private:
   };

   template <typename TSimType> ETDIteration<TSimType>::ETDIteration()
   {
   }

}

#endif // ETDITERATION_HPP
