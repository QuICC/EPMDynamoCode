/** \file SetProduct.hpp
 */

#ifndef SETPRODUCT_HPP
#define SETPRODUCT_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "General/EPMTypedefs.hpp"

namespace EPMDynamo {

   /**
    * @brief Set the results of the matrix multiplication
    *
    * \bug Needs cleaning and restructuring
    */
   class SetProduct
   {
      public:
         /**
          * @brief Set result
          */
         template<typename T1, typename T2>  static void eval(T1 &rOutVal, const T2 &inVal);

         /**
          * @brief Set result
          */
         template<typename T>  static void eval(MatrixZ &rOutVal, const int k, const T &inVal);
         
      protected:

      private:
   };

   template<typename T1, typename T2> inline void SetProduct::eval(T1 &rOutVal, const T2 &inVal)
   {
      rOutVal = inVal;
   }

   template<typename T> inline void SetProduct::eval(MatrixZ &rOutVal, const int k, const T &inVal)
   {
      rOutVal.col(k) = inVal;
   }

}

#endif // SETPRODUCT_HPP
