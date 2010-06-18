/** \file ComputedOperator.hpp
 */

#ifndef COMPUTEDOPERATOR_HPP
#define COMPUTEDOPERATOR_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "General/EPMTypedefs.hpp"

namespace EPMDynamo {

   /**
    * @brief This class implements the generic polynomial operator
    *
    * \bug Needs cleaning and restructuring
    */
   template <typename TExp> class ComputedOperator : public TExp
   {
      public:
         /**
         * @brief Constructor for a P operator
         */
         ComputedOperator(const Matrix &poly);

         /**
         * @brief Constructor for a C P operator
         */
         ComputedOperator(const typename TExp::ConstantType c, const Matrix &poly);

         /**
         * @brief Constructor for a RP or LP operator
         */
         ComputedOperator(const Array &arr, const Matrix &poly);

         /**
         * @brief Constructor for a CRP or CLP operator
         */
         ComputedOperator(const typename TExp::ConstantType c, const Array &arr, const Matrix &poly);

         /**
         * @brief Constructor for a LRP operator
         */
         ComputedOperator(const Array &left, const Array &right, const Matrix &poly);

         /**
         * @brief Constructor for a CLRP operator
         */
         ComputedOperator(const typename TExp::ConstantType c, const Array &left, const Array &right, const Matrix &poly);

         /**
         * @brief Destructor
         */
         virtual ~ComputedOperator();
         
      protected:
         /**
         * @brief Constructor
         */
         ComputedOperator();

      private:

   };

   template <typename TExp> ComputedOperator<TExp>::ComputedOperator(const Matrix &poly)
      : TExp(poly)
   {
   }

   template <typename TExp> ComputedOperator<TExp>::ComputedOperator(const typename TExp::ConstantType c, const Matrix &poly)
      : TExp(c, poly)
   {
   }

   template <typename TExp> ComputedOperator<TExp>::ComputedOperator(const Array &arr, const Matrix &poly)
      : TExp(arr, poly)
   {
   }

   template <typename TExp> ComputedOperator<TExp>::ComputedOperator(const typename TExp::ConstantType c, const Array &arr, const Matrix &poly)
      : TExp(c, arr, poly)
   {
   }

   template <typename TExp> ComputedOperator<TExp>::ComputedOperator(const Array &left, const Array &right, const Matrix &poly)
      : TExp(left, right, poly)
   {
   }

   template <typename TExp> ComputedOperator<TExp>::ComputedOperator(const typename TExp::ConstantType c, const Array &left, const Array &right, const Matrix &poly)
      : TExp(c, left, right, poly)
   {
   }

   template <typename TExp> ComputedOperator<TExp>::~ComputedOperator()
   {
   }

}

#endif // COMPUTEDOPERATOR_HPP
