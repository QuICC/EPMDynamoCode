/** \file ETDNOperators.hpp
 *  \brief Implementation of the ETDN operators
 */

#ifndef ETDNOPERATORS_HPP
#define ETDNOPERATORS_HPP

// System includes
//
#include <vector>

// External includes
//

// Project includes
//
#include "Domain/Truncation.hpp"
#include "Operators/DenseOperator.hpp"
#include "Timestepping/TimestepParameters.hpp"
#include "Timestepping/ETD/ETDOperators.hpp"
#include "Timestepping/ETD/ETDSchemeTraits.hpp"
#include "Simulations/Traits/SimulationTraits.hpp"
#include "LAPACK_Iface.hpp"

namespace EPMDynamo {

   /**
    * \brief Implementation of the ETDN operators
    */
   template <typename TSimType, int TSchemeOrder> class ETDNOperators
   {
      public:
         /// Typedef from Simulation trait to local radial basis type
         typedef typename TSimType::RadialBasisType    BasisType;

         /// Typedef from Simulation trait to local scalar type
         typedef typename TSimType::ScalarType    ScalarType;

         /// Typedef for an ETD Operator
         typedef typename ETDSchemeTraits<TSimType>::Operators ETDOps;

         /// Typedef for a smart pointer to ETDOperators
         typedef EPMSHARED_PTR<ETDOps> SmartETDOperators;

         /**
          * @brief Constructor
          *
          * @param pTrunc Truncation information
          */
         ETDNOperators(SmartTruncation pTrunc);

         /**
          * @brief Destructor
          */
         virtual ~ETDNOperators() {};

         /**
          * @brief Update the timestep matrices
          *
          * Still a pure virtual function
          *
          * @param dt   New timestep value
          * @param basis Radial basis
          */
         virtual void update(const EPMFloat dt, const BasisType &basis) = 0;

         /**
          * @brief Init the operators
          */
         void initOperators();

         /**
          * @brief Add a boundary condition
          *
          * @param pBC Boundary condition
          */
         void addBC(SmartBC pBC);

         /**
          * @brief Update the ETDN operators
          *
          * @param c Maximum eigen value
          * @param basis Radial basis
          */
         virtual void createOperators(const EPMFloat c, const BasisType &basis) = 0;

         /**
          * @brief Get the Fn operators
          *
          * @param n Index of the operator
          */
         const ETDOps& etdF(const int n) const;

         /**
          * @brief Get smart pointer to Fn operators
          *
          * @param n Index of the operator
          */
         SmartETDOperators pEtdF(const int n) const;

         /**
          * @brief Set the Fn operators
          *
          * @param n Index of the operator
          */
         ETDOps& rEtdF(const int n);


      protected:
         /**
          * @brief Update the number of required scalings
          *
          * @param c New number of scalings
          */
         void updateScalings(const EPMFloat c);

         /**
          * @brief Compute the scaled F0 values
          */
         void computeScaledF0();

         /**
          * @brief Compute the squared F0 values
          */
         void computeSquaredF0();

         /**
          * @brief Compute the squared F0 and F1 values
          */
         void computeSquaredF1();

         /**
          * @brief Compute the squared F0, F1 and F2 values
          */
         void computeSquaredF2();

         /**
          * @brief Compute the squared F0, F1, F2 and F3 values
          */
         void computeSquaredF3();

         /**
          * @brief Compute of the given matrix
          *
          * @param rMat Matrix to invert
          */
         void computeInverse(Matrix &rMat);

      private:
         /**
          * @brief Threshold value for requiring the scaling/squaring approach
          */
         static const EPMFloat  SCALINGSQUARING_THRESHOLD;

         /**
          * @brief Number of operators
          */
         const int   mNOps;

         /**
          * @brief Number of scalings/squaring required
          */
         int   mScalings;

         /**
          * @brief Absolute value of maximum eigenvalue
          */
         EPMFloat mMaxEig;

         /**
          * @brief Smart truncation information
          */
         SmartTruncation mpTrunc;

         /**
          * @brief Storage for the \f$M_i\f$ matrices required in the ETDN schemes
          */
         std::vector<SmartETDOperators>  mOperators;

         /**
          * @brief Compute the exponential of the given matrix
          *
          * @param mat Matrix to exponentiate
          *
          * \epmBug Currently uses a naive Taylor series
          */
         void computeExponential(Matrix& mat);

         /**
          * @brief Compute the F0(2z)
          */
         void squareF0();

         /**
          * @brief Compute the F1(2z)
          */
         void squareF1();

         /**
          * @brief Compute the F2(2z)
          */
         void squareF2();

         /**
          * @brief Compute the F3(2z)
          */
         void squareF3();

         /**
          * @brief Initialise the vector of ETDOperators
          */
         void initStorage();
   };

   template <typename TSimType, int TSchemeOrder> inline const typename ETDNOperators<TSimType, TSchemeOrder>::ETDOps& ETDNOperators<TSimType, TSchemeOrder>::etdF(const int n) const
   {
      return *(this->mOperators.at(n));
   }

   template <typename TSimType, int TSchemeOrder> inline  typename ETDNOperators<TSimType, TSchemeOrder>::SmartETDOperators ETDNOperators<TSimType, TSchemeOrder>::pEtdF(const int n) const
   {
      return this->mOperators.at(n);
   }

   template <typename TSimType, int TSchemeOrder> inline  typename ETDNOperators<TSimType, TSchemeOrder>::ETDOps& ETDNOperators<TSimType, TSchemeOrder>::rEtdF(const int n)
   {
      return *(this->mOperators.at(n));
   }

   template <typename TSimType, int TSchemeOrder> ETDNOperators<TSimType, TSchemeOrder>::ETDNOperators(SmartTruncation pTrunc)
      : mNOps(TSchemeOrder), mScalings(0), mMaxEig(0), mpTrunc(pTrunc)
   {
      // Initialise operators
      this->initStorage();
   }

   template <typename TSimType, int TSchemeOrder> void ETDNOperators<TSimType, TSchemeOrder>::addBC(SmartBC pBC)
   {
      // add boundary condition to all operators
      for(int i=0; i < this->mNOps; ++i)
      {
         this->mOperators.at(i)->addBC(pBC);
      }
   }

   template <typename TSimType, int TSchemeOrder> void ETDNOperators<TSimType, TSchemeOrder>::initOperators()
   {
      // init all operators
      for(int i=0; i < this->mNOps; ++i)
      {
         this->mOperators.at(i)->initOperators();
      }
   }

   template <typename TSimType, int TSchemeOrder> void ETDNOperators<TSimType, TSchemeOrder>::initStorage()
   {
      SmartETDOperators pOp;

      for(int i=0; i < this->mNOps; ++i)
      {
         pOp = SmartETDOperators(new ETDOps(this->mpTrunc));

         this->mOperators.push_back(pOp);
      }
   }

   template <typename TSimType, int TSchemeOrder> const EPMFloat  ETDNOperators<TSimType, TSchemeOrder>::SCALINGSQUARING_THRESHOLD = 10;

   template <typename TSimType, int TSchemeOrder> void ETDNOperators<TSimType, TSchemeOrder>::updateScalings(const EPMFloat c)
   {
      this->mScalings = std::ceil(std::log(SCALINGSQUARING_THRESHOLD*c*this->mMaxEig)/std::log(2.0));
   }

   template <typename TSimType, int TSchemeOrder> void ETDNOperators<TSimType, TSchemeOrder>::computeInverse(Matrix &rMat)
   {
      int rows = rMat.rows();
      int cols = rMat.cols();

      int info;
      ArrayI ipiv(rows);

      // Call LAPACK dgetrf routine for factorisation
      dgetrf_(&rows, &cols, rMat.data(), &rows, ipiv.data(), &info);

      // Test success of computation through assert
      assert(info == 0);

      int lwork = -1;
      Array work(rows);

      // Call LAPACK dgetri routine for optimal work size
      dgetri_(&rows, rMat.data(), &rows, ipiv.data(), work.data(), &lwork, &info);

      // Test success of computation through assert
      assert(info == 0);

      // Set to optimal work size
      lwork = static_cast<int>(work(0));
      work.resize(lwork);

      // Call LAPACK dgetri routine for inverse computation
      dgetri_(&rows, rMat.data(), &rows, ipiv.data(), work.data(), &lwork, &info);

      // Test success of computation through assert
      assert(info == 0);
   }

   template <typename TSimType, int TSchemeOrder> void ETDNOperators<TSimType, TSchemeOrder>::computeScaledF0()
   {
      // Loop over all degrees
      for(int i = 0; i < this->etdF(0).nOp(); ++i)
      {
         // Rescale operator
         this->rEtdF(0).rHarmOp(i).rOp() *= std::pow(2.0, this->mScalings);

         // Compute exponential of operator
         this->computeExponential(this->rEtdF(0).rHarmOp(i).rOp());
      }
   }

   template <typename TSimType, int TSchemeOrder> void ETDNOperators<TSimType, TSchemeOrder>::computeExponential(Matrix& rMat)
   {
      // Add identity
      rMat.diagonal().cwise() += 1.0;

      // Storage for the factorial factor
      EPMFloat factor = 1.0;

      // Loop over a certain number of expansion factors
      for(int i=2; i < 12; ++i)
      {
         factor *= static_cast<EPMFloat>(i);

         rMat += (rMat*rMat)/factor;
      }
   }

   template <typename TSimType, int TSchemeOrder> void ETDNOperators<TSimType, TSchemeOrder>::squareF0()
   {
      // Loop over all degrees
      for(int i = 0; i < this->etdF(0).nOp(); ++i)
      {
         // Define homogeneous operator
         this->rEtdF(0).rHarmOp(i).rOp() *= this->etdF(0).harmOp(i).op();
      }
   }

   template <typename TSimType, int TSchemeOrder> void ETDNOperators<TSimType, TSchemeOrder>::computeSquaredF0()
   {
      for(int i = 0; i < this->mScalings; ++i)
      {
         // Square F0
         this->squareF0();
      }
   }

   template <typename TSimType, int TSchemeOrder> void ETDNOperators<TSimType, TSchemeOrder>::squareF1()
   {
      // Loop over all degrees
      for(int i = 0; i < this->etdF(0).nOp(); ++i)
      {
         this->rEtdF(1).rHarmOp(i).rOp() += this->etdF(1).harmOp(i).op()*this->etdF(0).harmOp(i).op();
         this->rEtdF(1).rHarmOp(i).rOp() *= 0.5;
      }
   }

   template <typename TSimType, int TSchemeOrder> void ETDNOperators<TSimType, TSchemeOrder>::computeSquaredF1()
   {
      for(int i = 0; i < this->mScalings; ++i)
      {
         // "Square" F1
         this->squareF1();

         // Square F0
         this->squareF0();
      }
   }

   template <typename TSimType, int TSchemeOrder> void ETDNOperators<TSimType, TSchemeOrder>::squareF2()
   {
      // Loop over all degrees
      for(int i = 0; i < this->etdF(0).nOp(); ++i)
      {
         // Define homogeneous operator
         this->rEtdF(2).rHarmOp(i).rOp() *= 2.0;
         this->rEtdF(2).rHarmOp(i).rOp() += this->etdF(1).harmOp(i).op()*this->etdF(1).harmOp(i).op();
         this->rEtdF(2).rHarmOp(i).rOp() *= 0.25;
      }
   }

   template <typename TSimType, int TSchemeOrder> void ETDNOperators<TSimType, TSchemeOrder>::computeSquaredF2()
   {
      for(int i = 0; i < this->mScalings; ++i)
      {
         // "Square" F2
         this->squareF2();

         // "Square" F1
         this->squareF1();

         // Square F0
         this->squareF0();
      }
   }

   template <typename TSimType, int TSchemeOrder> void ETDNOperators<TSimType, TSchemeOrder>::squareF3()
   {
      // Loop over all degrees
      for(int i = 0; i < this->etdF(0).nOp(); ++i)
      {
         // Define homogeneous operator
         this->rEtdF(3).rHarmOp(i).rOp().cwise() *= 2.0;
         this->rEtdF(3).rHarmOp(i).rOp() += this->etdF(1).harmOp(i).op()*this->etdF(2).harmOp(i).op() + this->etdF(2).harmOp(i).op();
         this->rEtdF(3).rHarmOp(i).rOp().cwise() *= 0.125;
      }
   }

   template <typename TSimType, int TSchemeOrder> void ETDNOperators<TSimType, TSchemeOrder>::computeSquaredF3()
   {
      for(int i = 0; i < this->mScalings; ++i)
      {
         // "Square" F3
         this->squareF3();

         // "Square" F2
         this->squareF2();

         // "Square" F1
         this->squareF1();

         // Square F0
         this->squareF0();
      }
   }

}

#endif // ETDNOPERATORS_HPP
