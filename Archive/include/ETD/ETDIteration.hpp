/** \file ETDIteration.hpp
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

   template <typename TSim> class ETDIteration
   {
      public:
         /// Typedef from Simulation trait to local scalar type
         typedef typename TSim::ScalarType    ScalarType;

         /**
          * @brief Constructor
          */
         ETDIteration();

         /**
          * @brief Destructor
          */
         virtual ~ETDIteration();

         /**
          * @brief Compute the implemented iteration
          */
         virtual void compute(ScalarType &rVar, ScalarType rNTerms) = 0;
         
      protected:

      private:
   };

   template <typename TSim> ETDIteration<TSim>::ETDIteration()
   {
   }

   template <typename TSim> ETDIteration<TSim>::~ETDIteration()
   {
   }

}

#endif // ETDITERATION_HPP
