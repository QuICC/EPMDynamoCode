/** \file SubProduct.hpp
 */

#ifndef SUBPRODUCT_HPP
#define SUBPRODUCT_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "General/EPMTypedefs.hpp"

namespace EPMDynamo {

   /**
    * @brief Substract the results of the matrix multiplication
    *
    * \epmBug Needs cleaning and restructuring
    */
   class SubProduct
   {
      public:
         /**
          * @brief Substract result
          */
         template<typename T1, typename T2>  static void eval(T1 &rOutVal, const T2 &inVal);

         /**
          * @brief Substract result
          */
         template<typename T>  static void eval(MatrixZ &rOutVal, const int k, const T &inVal);

      protected:

      private:
   };

   template<typename T1, typename T2> inline void SubProduct::eval(T1 &rOutVal, const T2 &inVal)
   {
      rOutVal -= inVal;
   }

   template<typename T> inline void SubProduct::eval(MatrixZ &rOutVal, const int k, const T &inVal)
   {
      rOutVal.col(k) -= inVal;
   }

}

#endif // SUBPRODUCT_HPP
