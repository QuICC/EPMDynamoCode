/** \file SchemeIteration.hpp
 *  \brief General implementation of an iterative scheme iteration (step)
 */

#ifndef SCHEMEITERATION_HPP
#define SCHEMEITERATION_HPP

// System includes
//

// External includes
//

// Project includes
//

namespace EPMDynamo {

   /**
    * \brief General implementation of an iterative scheme iteration (step)
    *
    * \tparam TSimType Type of the simulation
    */
   template <typename TSimType> class SchemeIteration
   {
      public:
         /// Typedef from Simulation trait to local scalar type
         typedef typename TSimType::ScalarType    ScalarType;

         /**
          * @brief Constructor
          */
         SchemeIteration();

         /**
          * @brief Destructor
          */
         virtual ~SchemeIteration() {};

         /**
          * @brief Compute the implemented iteration
          *
          * @param rVar Input/Output variable
          * @param rNTerms New non linear terms
          */
         virtual void compute(ScalarType &rVar, ScalarType &rNTerms) = 0;

         /**
          * @brief Does this iteration provide timestep error information
          */
         bool providesError() const;
         
      protected:

      private:
   };

   template <typename TSimType> SchemeIteration<TSimType>::SchemeIteration()
   {
   }

   template <typename TSimType> bool SchemeIteration<TSimType>::providesError() const
   {
      return false;
   }
}

#endif // SCHEMEITERATION_HPP
