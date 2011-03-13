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
          * @param c Stiffness constant
          * @param pTrunc Truncation information
          * @param hasL0 Is l=0 mode required?
          */
         ETDNOperators(EPMFloat c, SmartTruncation pTrunc, bool hasL0);

         /**
          * @brief Destructor
          */
         virtual ~ETDNOperators() {};

         /**
          * @brief Update the timestep matrices
          *
          * Still a pure virtual function
          *
          * @param dt New timestep value
          * @param basis Radial basis
          */
         virtual void update(const EPMFloat dt, const BasisType &basis) = 0;

         /**
          * @brief Init the operators
          */
         void initOperators(const BasisType &basis);

         /**
          * @brief Add a boundary condition
          *
          * @param pBC Boundary condition
          */
         void addBC(SmartBC pBC);

         /**
          * @brief Update the ETDN operators
          *
          * @param h Timestep length
          * @param basis Radial basis
          */
         virtual void createOperators(const EPMFloat h, const BasisType &basis) = 0;

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
          * @brief get stiffness constant c
          */
         EPMFloat c() const;

         /**
          * @brief Check if matrix is full rank
          */
         bool isFullRank(const int l) const;

         /**
          * @brief Update the number of required scalings
          *
          * @param h Timestep lenght
          */
         void updateScalings(const EPMFloat h);

         /**
          * @brief Scale given operator
          */
         void scaleOperator(Matrix& rMat, const int l);

         /**
          * @brief Compute the scaled F0 values
          */
         void computeScaledF0(const int l);

         /**
          * @brief Compute the scaled F0 and F1 values
          */
         void computeScaledF1(const int l);

         /**
          * @brief Compute the scaled F0, F1 and F2 values
          */
         void computeScaledF2(const int l);

         /**
          * @brief Compute the scaled F0, F1, F2 and F3 values
          */
         void computeScaledF3(const int l);

         /**
          * @brief Compute the squared F0 values
          */
         void computeSquaredF0(const int l);

         /**
          * @brief Compute the squared F0 and F1 values
          */
         void computeSquaredF1(const int l);

         /**
          * @brief Compute the squared F0, F1 and F2 values
          */
         void computeSquaredF2(const int l);

         /**
          * @brief Compute the squared F0, F1, F2 and F3 values
          */
         void computeSquaredF3(const int l);

         /**
          * @brief Compute the inverse of the given matrix
          *
          * @param rMat Matrix to invert
          */
         void computeInverse(Matrix &rMat, const int l);

      private:
         /**
          * @brief Threshold value for requiring the scaling/squaring approach
          */
         static const EPMFloat  SCALINGSQUARING_THRESHOLD;

         /**
          * @brief Storage for the stiffness constant c
          */
         EPMFloat mC;

         /**
          * @brief Is l=0 mode required?
          */
         bool mHasL0;

         /**
          * @brief Number of operators
          */
         const int   mNOps;

         /**
          * @brief Number of scalings/squaring required
          */
         ArrayI   mScalings;

         /**
          * @brief Is the matrix full rank?
          */
         ArrayB mFullRank;

         /**
          * @brief Absolute value of maximum eigenvalues of the whole set of operators
          */
         Array mMaxEig;

         /**
          * @brief Smart truncation information
          */
         SmartTruncation mpTrunc;

         /**
          * @brief Storage for the \f$M_i\f$ matrices required in the ETDN schemes
          */
         std::vector<SmartETDOperators>  mOperators;

         /**
          * @brief Compute the taylor series of the Fk functions
          *
          * @param mat Matrix to expand
          */
         void computeFkTaylor(Matrix& mat, int k);

         /**
          * @brief Compute the F0(2z)
          *
          * @param l Harmonic degree
          */
         void squareF0(const int l);

         /**
          * @brief Compute the F1(2z)
          */
         void squareF1(const int l);

         /**
          * @brief Compute the F2(2z)
          *
          * @param l Harmonic degree
          */
         void squareF2(const int l);

         /**
          * @brief Compute the F3(2z)
          *
          * @param l Harmonic degree
          */
         void squareF3(const int l);

         /**
          * @brief Initialise the vector of ETDOperators
          */
         void initStorage();

         /**
          * @brief Compute and set the maximum eigenvalues
          *
          * @param basis Radial basis
          */
         void initEigenvalues(const BasisType &basis);

         /**
          * @brief Compute the rank of the operators
          *
          * @param basis Radial basis
          */
         void initRanks(const BasisType &basis);
   };

   template <typename TSimType, int TSchemeOrder> inline EPMFloat ETDNOperators<TSimType, TSchemeOrder>::c() const
   {
      return this->mC;
   }

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

   template <typename TSimType, int TSchemeOrder> inline bool ETDNOperators<TSimType, TSchemeOrder>::isFullRank(const int l) const
   {
      return this->mFullRank(l);
   }

   template <typename TSimType, int TSchemeOrder> ETDNOperators<TSimType, TSchemeOrder>::ETDNOperators(EPMFloat c, SmartTruncation pTrunc, bool hasL0)
      : mC(c), mHasL0(hasL0), mNOps(TSchemeOrder), mScalings(pTrunc->local()->spec()->nL()), mFullRank(pTrunc->local()->spec()->nL()), mMaxEig(pTrunc->local()->spec()->nL()), mpTrunc(pTrunc)
   {
      // Initialise the scalings and eigen values
      this->mScalings.setConstant(0);
      this->mFullRank.setConstant(true);
      this->mMaxEig.setConstant(0.0);

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

   template <typename TSimType, int TSchemeOrder> void ETDNOperators<TSimType, TSchemeOrder>::initOperators(const ETDNOperators<TSimType, TSchemeOrder>::BasisType &basis)
   {
      // init all operators
      for(int i=0; i < this->mNOps; ++i)
      {
         this->mOperators.at(i)->initOperators();
      }

      // Compute the eigenvalues
      this->initRanks(basis);

      // Compute the eigenvalues
      this->initEigenvalues(basis);

   }

   template <typename TSimType, int TSchemeOrder> void ETDNOperators<TSimType, TSchemeOrder>::initRanks(const ETDNOperators<TSimType, TSchemeOrder>::BasisType &basis)
   {
      // Loop over all degrees
      for(int i = this->etdF(0).minL(); i < this->etdF(0).nOp(); ++i)
      {
         // Define homogeneous operator
         this->rEtdF(0).rHarmOp(i).constructBOperator(1.0, basis.at(i).specLaplacian());

         char jobu = 'N';
         char jobvt = 'N';

         int N = this->rEtdF(0).rHarmOp(i).nTau();
         int lwork = 5*N;
         int info;

         Array s(N);
         Array work(lwork);

         // Call LAPACK dgetrf routine for factorisation
         dgesvd_(&jobu, &jobvt, &N, &N, this->rEtdF(0).rHarmOp(i).rOp().data(), &N, s.data(), NULL, &N, NULL, &N, work.data(), &lwork, &info);

         // Test success of computation through assert
         assert(info == 0);

         // Set the rank flag of the operator, this information will be used in a further step
         if(s.array().minCoeff() == 0.0)
         {
            this->mFullRank(i) = false;
         }
      }
   }

   template <typename TSimType, int TSchemeOrder> void ETDNOperators<TSimType, TSchemeOrder>::initEigenvalues(const ETDNOperators<TSimType, TSchemeOrder>::BasisType &basis)
   {
      // Loop over all degrees
      for(int i = this->etdF(0).minL(); i < this->etdF(0).nOp(); ++i)
      {
         // Define homogeneous operator
         this->rEtdF(0).rHarmOp(i).constructBOperator(this->c(), basis.at(i).specLaplacian());

         char jobvl = 'N';
         char jobvr = 'N';

         int N = this->rEtdF(0).rHarmOp(i).nTau();
         int lwork = 3*N;
         int info;

         Array wr(N);
         Array wi(N);
         Array work(lwork);

         // Call LAPACK dgetrf routine for factorisation
         dgeev_(&jobvl, &jobvr, &N, this->rEtdF(0).rHarmOp(i).rOp().data(), &N, wr.data(), wi.data(), NULL, &N, NULL, &N, work.data(), &lwork, &info);

         // Test success of computation through assert
         assert(info == 0);

         // Test that all eigenvalues are real!
         assert(wi.sum() == 0.0);

         // Test that all eigenvalues are negative else you risk a nasty blowup!
         assert(wr.maxCoeff() <= 0.0);

         // Set the maximum eigenvalue
         this->mMaxEig(i) = wr.array().abs().maxCoeff();
      }
   }

   template <typename TSimType, int TSchemeOrder> void ETDNOperators<TSimType, TSchemeOrder>::initStorage()
   {
      SmartETDOperators pOp;

      for(int i=0; i < this->mNOps; ++i)
      {
         pOp = SmartETDOperators(new ETDOps(this->c(), this->mpTrunc, this->mHasL0));

         this->mOperators.push_back(pOp);
      }
   }

   template <typename TSimType, int TSchemeOrder> const EPMFloat  ETDNOperators<TSimType, TSchemeOrder>::SCALINGSQUARING_THRESHOLD = 10;

   template <typename TSimType, int TSchemeOrder> void ETDNOperators<TSimType, TSchemeOrder>::updateScalings(const EPMFloat h)
   {
      // Loop over all degrees
      for(int i = this->etdF(0).minL(); i < this->etdF(0).nOp(); ++i)
      {
         this->mScalings(i) = std::max(0.0, 1+std::ceil(std::log(SCALINGSQUARING_THRESHOLD*h*this->mMaxEig(i))/std::log(2.0)));
      }
      std::cerr << this->mScalings.transpose() << std::endl;

      // Loop over all the ETD operators
      for(int i=0; i < this->mNOps; ++i)
      {
         this->mOperators.at(i)->updateTimestep(h);
      }
   }

   template <typename TSimType, int TSchemeOrder> void ETDNOperators<TSimType, TSchemeOrder>::computeInverse(Matrix &rMat, const int l)
   {
      // Comput the inverse
      if(this->isFullRank(l))
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

      // Compute pseudoinverse
      } else
      {
         // Define homogeneous operator
         char jobu = 'A';
         char jobvt = 'A';

         int N = rMat.cols();
         int lwork = -1;
         int info;

         Array s(N);
         Array work(N);

         Matrix U(N,N);
         Matrix VT(N,N);

         // Call LAPACK dgesvd with lwork=-1 to get optimal size
         dgesvd_(&jobu, &jobvt, &N, &N, rMat.data(), &N, s.data(), U.data(), &N, VT.data(), &N, work.data(), &lwork, &info);

         // Test success of computation through assert
         assert(info == 0);

         // Set to optimal work size
         lwork = static_cast<int>(work(0));
         work.resize(lwork);

         // Call LAPACK dgesvd routine for SVD decomposition
         dgesvd_(&jobu, &jobvt, &N, &N, rMat.data(), &N, s.data(), U.data(), &N, VT.data(), &N, work.data(), &lwork, &info);

         // Test success of computation through assert
         assert(info == 0);

         // Create array of reciprocal values
         for(int i=0;i< N; i++)
         {
            if(s(i) != 0.0)
            {
               s(i) = 1.0/s(i);
            }
         }

         // Compute the pseudo inverse
         rMat = VT.transpose() * s.asDiagonal() * U.transpose();
      }
   }

   template <typename TSimType, int TSchemeOrder> void ETDNOperators<TSimType, TSchemeOrder>::scaleOperator(Matrix& rMat, const int l)
   {
         // Rescale operator
         rMat *= std::pow(2.0, -this->mScalings(l));
   }

   template <typename TSimType, int TSchemeOrder> void ETDNOperators<TSimType, TSchemeOrder>::computeScaledF0(const int l)
   {
      // Rescale operator
      this->scaleOperator(this->rEtdF(0).rHarmOp(l).rOp(), l);

      // Compute taylor series of F0
      this->computeFkTaylor(this->rEtdF(0).rHarmOp(l).rOp(), 0);
   }

   template <typename TSimType, int TSchemeOrder> void ETDNOperators<TSimType, TSchemeOrder>::computeScaledF1(const int l)
   {
      // Rescale operator
      this->scaleOperator(this->rEtdF(0).rHarmOp(l).rOp(), l);

      // Copy rescaled operator in F1
      this->rEtdF(1).rHarmOp(l).rOp() = this->etdF(0).harmOp(l).op();

      // Compute taylor series of F0
      this->computeFkTaylor(this->rEtdF(0).rHarmOp(l).rOp(), 0);

      // Compute taylor series of F1
      this->computeFkTaylor(this->rEtdF(1).rHarmOp(l).rOp(), 1);
   }

   template <typename TSimType, int TSchemeOrder> void ETDNOperators<TSimType, TSchemeOrder>::computeScaledF2(const int l)
   {
      // Rescale operator
      this->scaleOperator(this->rEtdF(0).rHarmOp(l).rOp(), l);

      // Copy rescaled operator in F1
      this->rEtdF(1).rHarmOp(l).rOp() = this->etdF(0).harmOp(l).op();

      // Copy rescaled operator in F1
      this->rEtdF(2).rHarmOp(l).rOp() = this->etdF(0).harmOp(l).op();

      // Compute taylor series of F0
      this->computeFkTaylor(this->rEtdF(0).rHarmOp(l).rOp(), 0);

      // Compute taylor series of F1
      this->computeFkTaylor(this->rEtdF(1).rHarmOp(l).rOp(), 1);

      // Compute taylor series of F2
      this->computeFkTaylor(this->rEtdF(2).rHarmOp(l).rOp(), 2);
   }

   template <typename TSimType, int TSchemeOrder> void ETDNOperators<TSimType, TSchemeOrder>::computeScaledF3(const int l)
   {
      // Rescale operator
      this->scaleOperator(this->rEtdF(0).rHarmOp(l).rOp(), l);

      // Copy rescaled operator in F1
      this->rEtdF(1).rHarmOp(l).rOp() = this->etdF(0).harmOp(l).op();

      // Copy rescaled operator in F2
      this->rEtdF(2).rHarmOp(l).rOp() = this->etdF(0).harmOp(l).op();

      // Copy rescaled operator in F3
      this->rEtdF(3).rHarmOp(l).rOp() = this->etdF(0).harmOp(l).op();

      // Compute taylor series of F0
      this->computeFkTaylor(this->rEtdF(0).rHarmOp(l).rOp(), 0);

      // Compute taylor series of F1
      this->computeFkTaylor(this->rEtdF(1).rHarmOp(l).rOp(), 1);

      // Compute taylor series of F2
      this->computeFkTaylor(this->rEtdF(2).rHarmOp(l).rOp(), 2);

      // Compute taylor series of F3
      this->computeFkTaylor(this->rEtdF(3).rHarmOp(l).rOp(), 3);
   }

   template <typename TSimType, int TSchemeOrder> void ETDNOperators<TSimType, TSchemeOrder>::computeFkTaylor(Matrix& rMat, int k)
   {
      // Storage for the factorial factor
      EPMFloat factor = 1.0;

      // Compute the starting factorial factor
      for(int i=2; i <= k; i++)
      {
         factor *= static_cast<EPMFloat>(i);
      }

      // Store input matrix
      Matrix inOp = rMat;

      // Set starting values to identity matrix
      Matrix tmp(rMat.rows(), rMat.cols());
      tmp.setConstant(0.0);
      tmp.diagonal().setConstant(1.0);

      // Set staring value to Identity/factor
      rMat.setConstant(0.0);
      rMat.diagonal().setConstant(1.0/factor);

      // Do the taylor expansion loop
      for(int i=k+1; i < 16; ++i)
      {
         factor *= static_cast<EPMFloat>(i);

         tmp *= inOp;

         rMat += tmp/factor;
      }
   }

   template <typename TSimType, int TSchemeOrder> void ETDNOperators<TSimType, TSchemeOrder>::squareF0(const int l)
   {
      // Define homogeneous operator
      this->rEtdF(0).rHarmOp(l).rOp() *= this->etdF(0).harmOp(l).op();
   }

   template <typename TSimType, int TSchemeOrder> void ETDNOperators<TSimType, TSchemeOrder>::computeSquaredF0(const int l)
   {
      // Perform the squarings
      for(int i = 0; i < this->mScalings(l); ++i)
      {
         // Square F0
         this->squareF0(l);
      }
   }

   template <typename TSimType, int TSchemeOrder> void ETDNOperators<TSimType, TSchemeOrder>::squareF1(const int l)
   {
      this->rEtdF(1).rHarmOp(l).rOp() += this->etdF(0).harmOp(l).op()*this->etdF(1).harmOp(l).op();
      this->rEtdF(1).rHarmOp(l).rOp() *= 0.5;
   }

   template <typename TSimType, int TSchemeOrder> void ETDNOperators<TSimType, TSchemeOrder>::computeSquaredF1(const int l)
   {
      // Perform the squarings
      for(int i = 0; i < this->mScalings(l); ++i)
      {
         // "Square" F1
         this->squareF1(l);

         // Square F0
         this->squareF0(l);
      }
   }

   template <typename TSimType, int TSchemeOrder> void ETDNOperators<TSimType, TSchemeOrder>::squareF2(const int l)
   {
      // Define homogeneous operator
      this->rEtdF(2).rHarmOp(l).rOp() *= 2.0;
      this->rEtdF(2).rHarmOp(l).rOp() += this->etdF(1).harmOp(l).op()*this->etdF(1).harmOp(l).op();
      this->rEtdF(2).rHarmOp(l).rOp() *= 0.25;
   }

   template <typename TSimType, int TSchemeOrder> void ETDNOperators<TSimType, TSchemeOrder>::computeSquaredF2(const int l)
   {
      // Perform the squarings
      for(int i = 0; i < this->mScalings(l); ++i)
      {
         // "Square" F2
         this->squareF2(l);

         // "Square" F1
         this->squareF1(l);

         // Square F0
         this->squareF0(l);
      }
   }

   template <typename TSimType, int TSchemeOrder> void ETDNOperators<TSimType, TSchemeOrder>::squareF3(const int l)
   {
      // Define homogeneous operator
      this->rEtdF(3).rHarmOp(l).rOp().array() *= 2.0;
      this->rEtdF(3).rHarmOp(l).rOp() += this->etdF(1).harmOp(l).op()*this->etdF(2).harmOp(l).op() + this->etdF(2).harmOp(l).op();
      this->rEtdF(3).rHarmOp(l).rOp().array() *= 0.125;
   }

   template <typename TSimType, int TSchemeOrder> void ETDNOperators<TSimType, TSchemeOrder>::computeSquaredF3(const int l)
   {
      // Perform the squarings
      for(int i = 0; i < this->mScalings(l); ++i)
      {
         // "Square" F3
         this->squareF3(l);

         // "Square" F2
         this->squareF2(l);

         // "Square" F1
         this->squareF1(l);

         // Square F0
         this->squareF0(l);
      }
   }

}

#endif // ETDNOPERATORS_HPP
