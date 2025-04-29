/** \file ErrorComputer.hpp
 *  @brief Implementation of some tools used to compute errors in the tests
 */

#ifndef ERRORCOMPUTER_HPP
#define ERRORCOMPUTER_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "General/EPMTypedefs.hpp"

namespace EPMDynamo {

   /**
    * @brief Implementation of some tools used to compute errors in the tests
    *
    * \epmTSTodo Need to finish documentation
    */
   class ErrorComputer
   {
      public:
         /**
         * @brief Constructor
         */
         ErrorComputer();

         /**
         * @brief Destructor
         */
         virtual ~ErrorComputer() {};

         /**
          * @brief Maximum allowed error for passing check
          */
         static double MAX_ALLOWED_ERROR;

         /**
          * @brief This static double set the threshold for which to compute relative
          * errors
          */
         static double RELATIVE_MAXERROR;

         /**
          * @brief Get max error
          */
         double max() const;

         /**
          * @brief Get relative max error
          */
         double relMax() const;

         /**
          * @brief Reset stored errors
          */
         void resetErrors();

         /**
          * @brief Check error 
          */
         template<typename T, int R, int C> int checkPrecision(const Eigen::Matrix<T, R, C>  &val1, const Eigen::Matrix<T, R, C> &val2);

         /**
          * @brief Check error
          */
         template<typename CT, int R, int C> int checkPrecision(const Eigen::Matrix<std::complex<CT>, R, C>  &val1, const Eigen::Matrix<std::complex<CT>, R, C> &val2);
         
      protected:
         /**
          * @brief Compute the maximum of the error
          */
         template<typename T, int R, int C> double maxError(const Eigen::Matrix<T, R, C>  &val1, const Eigen::Matrix<T, R, C> &val2);

         /**
          * @brief Compute the maximum of the error for complex data
          */
         template<typename CT, int R, int C> double maxError(const Eigen::Matrix<std::complex<CT>, R, C>  &val1, const Eigen::Matrix<std::complex<CT>, R, C> &val2);

         /**
          * @brief Compute the maximum of the relative error
          */
         template<typename T, int R, int C> double maxRelError(const Eigen::Matrix<T, R, C>  &val1, const Eigen::Matrix<T, R, C> &val2);

         /**
          * @brief Compute the maximum of the relative error for complex data
          */
         template<typename CT, int R, int C> double maxRelError(const Eigen::Matrix<std::complex<CT>, R, C>  &val1, const Eigen::Matrix<std::complex<CT>, R, C> &val2);

      private:
         /**
          * @brief Storage for the maximum error
          */
         double mMax;

         /**
          * @brief Storage for the maximum relative error
          */
         double mRelMax;
   };

   inline double ErrorComputer::max() const
   {
      return this->mMax;
   }

   inline double ErrorComputer::relMax() const
   {
      return this->mRelMax;
   }

   template<typename T, int R, int C> int ErrorComputer::checkPrecision(const Eigen::Matrix<T, R, C>  &val1, const Eigen::Matrix<T, R, C> &val2)
   {
      double error = this->maxError(val1, val2);
      double relError = this->maxRelError(val1, val2);

      if(error > this->mMax)
      {
         this->mMax = error;
      }

      if(relError > this->mRelMax)
      {
         this->mRelMax = relError;
      }

      int status;

      if(relError > MAX_ALLOWED_ERROR)
      {
         status = 1;
      } else
      {
         status = 0;
      }

      return status;
   }

   template<typename CT, int R, int C> int ErrorComputer::checkPrecision(const Eigen::Matrix<std::complex<CT>, R, C>  &val1, const Eigen::Matrix<std::complex<CT>, R, C> &val2)
   {
      double error = this->maxError(val1, val2);
      double relError = this->maxRelError(val1, val2);

      if(error > this->mMax)
      {
         this->mMax = error;
      }

      if(relError > this->mRelMax)
      {
         this->mRelMax = relError;
      }

      int status;

      if(relError > MAX_ALLOWED_ERROR)
      {
         status = 1;
      } else
      {
         status = 0;
      }

      return status;
   }

   template<typename T, int R, int C> double ErrorComputer::maxError(const Eigen::Matrix<T, R, C>  &val1, const Eigen::Matrix<T, R, C> &val2)
   {
      double errMax = (val1-val2).array().abs().maxCoeff();

      return errMax;
   }

   template<typename CT, int R, int C> double ErrorComputer::maxError(const Eigen::Matrix<std::complex<CT>, R, C>  &val1, const Eigen::Matrix<std::complex<CT>, R, C> &val2)
   {
      Eigen::Matrix<std::complex<CT>, Eigen::Dynamic, Eigen::Dynamic> tmp(val1.rows(), val1.cols());

      tmp = val1 - val2;

      double errRe = tmp.real().array().abs().maxCoeff();
      double errIm = tmp.imag().array().abs().maxCoeff();

      double errMax = std::max(errRe, errIm);

      return errMax;
   }

   template<typename T, int R, int C> double ErrorComputer::maxRelError(const Eigen::Matrix<T, R, C>  &val1, const Eigen::Matrix<T, R, C> &val2)
   {
      Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> tmp(val1.rows(), val1.cols());

      tmp = val1 - val2;

      for(int j=0; j < val1.cols(); ++j)
      {
         for(int i=0; i < val1.rows(); ++i)
         {
            if(abs(val1(i,j)) > RELATIVE_MAXERROR)
            {
               tmp(i,j) = tmp(i,j)/val1(i,j);
            } else
            {
               tmp(i,j) = std::abs(tmp(i,j));
            }
         }
      }

      double errMax =  tmp.array().abs().maxCoeff();

      return errMax;
   }

   template<typename CT, int R, int C> double ErrorComputer::maxRelError(const Eigen::Matrix<std::complex<CT>, R, C>  &val1, const Eigen::Matrix<std::complex<CT>, R, C> &val2)
   {
      Eigen::Matrix<std::complex<CT>, Eigen::Dynamic, Eigen::Dynamic> tmp(val1.rows(), val1.cols());

      tmp = val1 - val2;

      for(int j=0; j < val1.cols(); ++j)
      {
         for(int i=0; i < val1.rows(); ++i)
         {
            if(std::abs(val1(i,j).real()) > RELATIVE_MAXERROR)
            {
               tmp(i,j).real() = tmp(i,j).real()/val1(i,j).real();
            } else{
               tmp(i,j).real() = std::abs(tmp(i,j).real());
            }
            if(std::abs(val1(i,j).imag()) > RELATIVE_MAXERROR)
            {
               tmp(i,j).imag() = tmp(i,j).imag()/val1(i,j).imag();
            } else{
               tmp(i,j).imag() = std::abs(tmp(i,j).imag());
            }
         }
      }

      double errRe = tmp.real().array().abs().maxCoeff();
      double errIm = tmp.imag().array().abs().maxCoeff();

      double errMax = std::max(errRe, errIm);

      return errMax;
   }

}

#endif // ERRORCOMPUTER_HPP
