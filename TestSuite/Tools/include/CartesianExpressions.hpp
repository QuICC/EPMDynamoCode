/** \file CartesianExpressions.hpp
 *  @brief Implementation of a simple way to express the fields in cartesian coordinates to simplfy production of test cases
 */

#ifndef CARTESIANEXPRESSIONS_HPP
#define CARTESIANEXPRESSIONS_HPP

// System includes
//
#include <Eigen/Array>

// Project includes
//
#include "General/EPMTypedefs.hpp"
#include "Domain/Truncation.hpp"

namespace EPMDynamo {

   /**
    * @brief Implementation of a simple way to express the fields in cartesian coordinates to simplfy production of test cases
    */
   class CartesianExpressions
   {
      public:
         /**
          * @brief Constructor
          *
          * @param pTrunc Truncation information
          */
         CartesianExpressions(SmartTruncation pTrunc);

         /**
          * @brief Simple empty destructor
          */
         virtual ~CartesianExpressions() {};

         /**
          * @brief X component
          *
          * @param n Radial index
          * @param power Power to compute
          */
         Matrix  x(const int n, const int power);

         /**
          * @brief Y component
          *
          * @param n Radial index
          * @param power Power to compute
          */
         Matrix  y(const int n, const int power);

         /**
          * @brief Z component
          *
          * @param n Radial index
          * @param power Power to compute
          */
         Matrix  z(const int n, const int power);

         /**
          * @brief XYZ to Radial component
          *
          * @param xComp X component
          * @param yComp Y component
          * @param zComp Z component
          * @param rIdx Radial index
          */
         Matrix xyzToR(const Matrix& xComp, const Matrix& yComp, const Matrix& zComp, const int rIdx);

         /**
          * @brief XYZ to Theta component
          *
          * @param xComp X component
          * @param yComp Y component
          * @param zComp Z component
          * @param rIdx Radial index
          */
         Matrix xyzToTh(const Matrix& xComp, const Matrix& yComp, const Matrix& zComp, const int rIdx);

         /**
          * @brief XYZ to Phi component
          *
          * @param xComp X component
          * @param yComp Y component
          * @param zComp Z component
          * @param rIdx Radial index
          */
         Matrix xyzToPh(const Matrix& xComp, const Matrix& yComp, const Matrix& zComp, const int rIdx);
         
      private:
         /**
          * @brief Store the spectral truncation object
          */
         SmartTruncation mpTrunc;

         /**
          * @brief Number phi points
          */
         int mNph;

         /**
          * @brief Store the radial grid
          */
         Array    mRad;

         /**
          * @brief Store the sin theta grid
          */
         Array    mSth;

         /**
          * @brief Store the cos theta grid
          */
         Array    mCth;

         /**
          * @brief Store the phi grid
          */
         Array    mPhi;

         /**
          * @brief Get the truncation information
          */
         SmartTruncation trunc() const;

         /**
          * @brief Get the X component matrix
          *
          * @param n Radial index
          */
         Matrix xMat(const int n);

         /**
          * @brief Get the Y component matrix
          *
          * @param n Radial index
          */
         Matrix yMat(const int n);

         /**
          * @brief Get the Z component matrix
          *
          * @param n Radial index
          */
         Matrix zMat(const int n);

         /**
          * @brief Get the X to R component conversion matrix
          *
          * @param n Radial index
          */
         Matrix xToRMat(const int n);

         /**
          * @brief Get the Y to R component conversion matrix
          *
          * @param n Radial index
          */
         Matrix yToRMat(const int n);

         /**
          * @brief Get the Z to R component conversion matrix
          *
          * @param n Radial index
          */
         Matrix zToRMat(const int n);

         /**
          * @brief Get the the X to Theta component conversion matrix
          *
          * @param n Radial index
          */
         Matrix xToThMat(const int n);

         /**
          * @brief Get the Y to Theta component conversion matrix
          *
          * @param n Radial index
          */
         Matrix yToThMat(const int n);

         /**
          * @brief Get the Z to Theta component conversion matrix
          *
          * @param n Radial index
          */
         Matrix zToThMat(const int n);

         /**
          * @brief Get the X to Phi component conversion matrix
          *
          * @param n Radial index
          */
         Matrix xToPhMat(const int n);

         /**
          * @brief Get the Y to Phi component conversion matrix
          *
          * @param n Radial index
          */
         Matrix yToPhMat(const int n);

         /**
          * @brief Get the Z to Phi component conversion matrix
          *
          * @param n Radial index
          */
         Matrix zToPhMat(const int n);

   };

   inline SmartTruncation CartesianExpressions::trunc() const
   {
      return this->mpTrunc;
   }

   inline Matrix  CartesianExpressions::x(const int n, const int power){
      Matrix tmp = std::pow(this->mRad(n), power)*this->xMat(n).cwise().pow(power);
      return tmp;
   };

   inline Matrix  CartesianExpressions::y(const int n, const int power){
      Matrix tmp = std::pow(this->mRad(n), power)*this->yMat(n).cwise().pow(power);
      return tmp;
   };

   inline Matrix  CartesianExpressions::z(const int n, const int power){
      Matrix tmp = std::pow(this->mRad(n), power)*this->zMat(n).cwise().pow(power);
      return tmp;
   };

   inline Matrix CartesianExpressions::xyzToR(const Matrix& xComp, const Matrix& yComp, const Matrix& zComp, const int rIdx)
   {
      Matrix tmp;
      tmp = xComp.cwise()*this->xToRMat(rIdx) + yComp.cwise()*this->yToRMat(rIdx) + zComp.cwise()*this->zToRMat(rIdx);

      return tmp;
   };

   inline Matrix CartesianExpressions::xyzToTh(const Matrix& xComp, const Matrix& yComp, const Matrix& zComp, const int rIdx)
   {
      Matrix tmp;
      tmp = xComp.cwise()*this->xToThMat(rIdx) + yComp.cwise()*this->yToThMat(rIdx) + zComp.cwise()*this->zToThMat(rIdx);

      return tmp;
   };

   inline Matrix CartesianExpressions::xyzToPh(const Matrix& xComp, const Matrix& yComp, const Matrix& zComp, const int rIdx)
   {
      Matrix tmp;
      tmp = xComp.cwise()*this->xToPhMat(rIdx) + yComp.cwise()*this->yToPhMat(rIdx) + zComp.cwise()*this->zToPhMat(rIdx);

      return tmp;
   };

}

#endif // CARTESIANEXPRESSIONS_HPP
