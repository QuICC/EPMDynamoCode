/** \file StoredOperator.hpp
 */

#ifndef STOREDOPERATOR_HPP
#define STOREDOPERATOR_HPP

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
   template <typename TExp> class StoredOperator : public TExp
   {
      public:
         /**
         * @brief Constructor for a RP or LP operator
         */
         StoredOperator(const Array &arr, const Matrix &poly);

         /**
         * @brief Constructor for a CRP or CLP operator
         */
         StoredOperator(const typename TExp::ConstantType c, const Array &arr, const Matrix &poly);

         /**
         * @brief Constructor for a LRP operator
         */
         StoredOperator(const Array &left, const Array &right, const Matrix &poly);

         /**
         * @brief Constructor for a CLRP operator
         */
         StoredOperator(const typename TExp::ConstantType c, const Array &left, const Array &right, const Matrix &poly);

         /**
         * @brief Destructor
         */
         virtual ~StoredOperator();

         /**
          * @brief Get stored product operator
          */
         const Matrix& productOp() const;

         /**
          * @brief Get operator
          */
         typename TExp::OutputMatrix op() const;
         
      protected:
         /**
         * @brief Constructor
         */
         StoredOperator();

      private:
         /**
          * @brief Storage for the operator
          */
         Matrix   mStoredOperator;

   };

   template <typename TExp> inline const Matrix& StoredOperator<TExp>::productOp() const
   {
      return mStoredOperator;
   }

   template <typename TExp> inline typename TExp::OutputMatrix StoredOperator<TExp>::op() const
   {
      return this->c()*this->productOp();
   }

   template <typename TExp> StoredOperator<TExp>::StoredOperator(const Array &arr, const Matrix &poly)
      : TExp(arr, poly)
   {
      mStoredOperator = TExp::productOp();
   }

   template <typename TExp> StoredOperator<TExp>::StoredOperator(const typename TExp::ConstantType c, const Array &arr, const Matrix &poly)
      : TExp(c, arr, poly)
   {
      mStoredOperator = TExp::productOp();
   }

   template <typename TExp> StoredOperator<TExp>::StoredOperator(const Array &left, const Array &right, const Matrix &poly)
      : TExp(left, right, poly)
   {
      mStoredOperator = TExp::productOp();
   }

   template <typename TExp> StoredOperator<TExp>::StoredOperator(const typename TExp::ConstantType c, const Array &left, const Array &right, const Matrix &poly)
      : TExp(c, left, right, poly)
   {
      mStoredOperator = TExp::productOp();
   }

   template <typename TExp> StoredOperator<TExp>::~StoredOperator()
   {
   }

}

#endif // STOREDOPERATOR_HPP
