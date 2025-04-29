/** \file AddProduct.hpp
 */

#ifndef ADDPRODUCT_HPP
#define ADDPRODUCT_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "General/EPMTypedefs.hpp"

namespace EPMDynamo {

   /**
    * @brief Add the results of the matrix multiplication
    */
   class AddProduct
   {
      public:
         /**
          * @brief Add result
          */
         template<typename T1, typename T2>  static void eval(T1 &rOutVal, const T2 &inVal);

         /**
          * @brief Add result
          */
         template<typename T>  static void eval(MatrixZ &rOutVal, const int k, const T &inVal);
      protected:

      private:
   };

   template<typename T1, typename T2> inline void AddProduct::eval(T1 &rOutVal, const T2 &inVal)
   {
      rOutVal += inVal;
   }

   template<typename T> inline void AddProduct::eval(MatrixZ &rOutVal, const int k, const T &inVal)
   {
      rOutVal.col(k) += inVal;
   }

}

#endif // ADDPRODUCT_HPP
