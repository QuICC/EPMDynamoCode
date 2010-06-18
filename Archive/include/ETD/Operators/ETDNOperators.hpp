/** \file ETDNOperators.hpp
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
#include "Domain/SpectralTruncation.hpp"
#include "Timestepping/TimestepParameters.hpp"
#include "Timestepping/ETD/ETDOperators.hpp"
#include "Simulations/SimulationTraits.hpp"

namespace EPMDynamo {

   template <typename TSim, int TNOps> class ETDNOperators
   {
      public:
         /// Typedef from Simulation trait to local radial basis type
         typedef typename TSim::RadialBasisType    BasisType;

         /// Typedef from Simulation trait to local scalar type
         typedef typename TSim::ScalarType    ScalarType;

         /// Typedef for a smart pointer to ETDOperators
         typedef typename ETDOperators<TSim> SmartETDOperators;

         /**
          * @brief Constructor
          *
          * @param a Coefficient \f$a\f$ of timestep scheme
          * @param b Coefficient \f$b\f$ of timestep scheme
          * @param basis Reference to the basis used for the operators
          */
         ETDNOperators(SmartSTrunc pTrunc);

         /**
          * @brief Destructor
          */
         virtual ~ETDNOperators();

         /**
          * @brief Update the ETDN operators
          */
         void createOperators(const DynamoFloat h, const BasisType &basis) = 0;

      protected:
         /**
          * @brief Get the Fn operators
          */
         const ETDOperators& etdF(const int n) const;

         /**
          * @brief Get smart pointer to Fn operators
          */
         SmartETDOperators pEtdF(const int n) const;

         /**
          * @brief Set the Fn operators
          */
         ETDOperators& rEtdF(const int n);

         /**
          * @brief Update the number of required scalings
          */
         void updateScalings(const DynamoFloat c);

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

      private:
         /**
          * @brief Threshold value for requiring the scaling/squaring approach
          */
         static const DynamoFloat  SCALINGSQUARING_THRESHOLD;

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
         DynamoFloat mMaxEig;

         /**
          * @brief Smart truncation information
          */
         SmartSTrunc mpSTrunc;

         /**
          * @brief Storage for the \f$M_i\f$ matrices required in the ETDN schemes
          */
         std::vector<SmartETDoperators>  mOperators;

         /**
          * @brief Compute the exponential of the given matrix
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
          * @brier Initialise the vector of ETDOperators
          */
         void initStorage();
   };

   template <typename TSim, int TNOps> inline const ETDOperators& ETDNOperators<TSim, TNOps>::etdF(const int n)
   {
      return *(this->mOperators.at(n));
   }

   template <typename TSim, int TNOps> inline SmartETDOperators ETDNOperators<TSim, TNOps>::pEtdF(const int n)
   {
      return this->mOperators.at(n);
   }

   template <typename TSim, int TNOps> inline ETDOperators& ETDNOperators<TSim, TNOps>::rEtdF(const int n)
   {
      return *(this->mOperators.at(n));
   }

   template <typename TSim, int TNOps> ETDNOperators<TSim, TNOps>::ETDNOperators(SmartSTrunc pTrunc)
      : mNOps(TNOps), mpSTrunc(pTrunc), mScalings(0), mMaxEig(0)
   {
      // Initialise operators
      this->initStorage();
   }

   template <typename TSim, int TNOps> ETDNOperators<TSim, TNOps>::~ETDNOperators()
   {
   }

   template <typename TSim, int TNOps> void ETDNOperators<TSim, TNOps>::initStorage()
   {
      SmartETDOperators pOp;

      for(int i=0; i < this->mNOps; ++i)
      {
         pOp = SmartETDOperators(new ETDOperators(this->mpSTrunc));

         this->mOperators.push_back(pOp);
      }
   }

   template <typename TSim, int TNOps> const DynamoFloat  ETDNOperators<TSim, TNOps>::SCALINGSQUARING_THRESHOLD = 10;

   template <typename TSim, int TNOps> void ETDNOperators<TSim, TNOps>::updateScalings(const DynamoFloat c)
   {
      this->mScalings = std::ceil(std::log(SCALINGSQUARING_THRESHOLD*c*this->mMaxEig)/std::log(2.0));
   }

   template <typename TSim, int TNOps> void ETDNOperators<TSim, TNOps>::computeScaledF0()
   {
      // Storage for "id" of operator
      int l;

      // Loop over all degrees
      for(int i = 0; i < this->etdF(0).nOp(); ++i)
      {
         // Get degree of the current operator
         l = this->rEtdF(0).op(i).id();

         // Rescale operator
         this->rEtdF(0).rOp(i).rOp().cwise() *= std::pow(2.0, this->mScalings);

         // Compute exponential of operator
         this->computeExponential(this->rEtdF(0).rOp(i).rOp());
      }
   }

   template <typename TSim, int TNOps> void ETDNOperators<TSim, TNOps>::computeExponential(Matrix& rMat)
   {
      // FOR THE MOMENT USE A BASIC TAYLOR SERIES

      // Add identity
      rMat.diagonal.cwise() += 1.0;

      // Storage for the factorial factor
      DynamoFloat factor = 1.0;

      // Loop over a certain number of expansion factors
      for(int i=2; i < 12; ++i)
      {
         factor *= static_cast<DynamoFloat>(i);

         rMat += rMat*rMat.cwise()/factor;
      }
   }

   template <typename TSim, int TNOps> void ETDNOperators<TSim, TNOps>::squareF0()
   {
      // Storage for "id" of operator
      int l;

      // Loop over all degrees
      for(int i = 0; i < this->etdF(0).nOp(); ++i)
      {
         // Get degree of the current operator
         l = this->rEtdF(0).op(i).id();

         // Define homogeneous operator
         this->rEtdF(0).rOp(i).rOp() *= this->etdF(0).op(i).op();
      }
   }

   template <typename TSim, int TNOps> void ETDNOperators<TSim, TNOps>::computeSsquaredF0()
   {
      for(int i = 0; i < this->mScalings; ++i)
      {
         // Square F0
         this->squareF0();
      }
   }

   template <typename TSim, int TNOps> void ETDNOperators<TSim, TNOps>::squareF1()
   {
      // Storage for "id" of operator
      int l;

      // Loop over all degrees
      for(int i = 0; i < this->etdF(0).nOp(); ++i)
      {
         // Get degree of the current operator
         l = this->rEtdF(0).op(i).id();

         this->rEtdF(1).rOp(i).rOp() += this->etdF(1).op(i).op()*this->etdF(0).op(i).op();
         this->rEtdF(1).rOp(i).rOp().cwise() *= 0.5;
      }
   }

   template <typename TSim, int TNOps> void ETDNOperators<TSim, TNOps>::computeSsquaredF1()
   {
      for(int i = 0; i < this->mScalings; ++i)
      {
         // "Square" F1
         this->squareF1();

         // Square F0
         this->squareF0();
      }
   }

   template <typename TSim, int TNOps> void ETDNOperators<TSim, TNOps>::squareF2()
   {
      // Storage for "id" of operator
      int l;

      // Loop over all degrees
      for(int i = 0; i < this->etdF(0).nOp(); ++i)
      {
         // Get degree of the current operator
         l = this->rEtdF(0).op(i).id();

         // Define homogeneous operator
         this->rEtdF(2).rOp(i).rOp().cwise() *= 2.0;
         this->rEtdF(2).rOp(i).rOp() += this->etdF(1).op(i).op()*this->etdF(1).op(i).op();
         this->rEtdF(2).rOp(i).rOp() *= 0.25;
      }
   }

   template <typename TSim, int TNOps> void ETDNOperators<TSim, TNOps>::computeSsquaredF2()
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

   template <typename TSim, int TNOps> void ETDNOperators<TSim, TNOps>::squareF3()
   {
      // Storage for "id" of operator
      int l;

      // Loop over all degrees
      for(int i = 0; i < this->etdF(0).nOp(); ++i)
      {
         // Get degree of the current operator
         l = this->rEtdF(0).op(i).id();

         // Define homogeneous operator
         this->rEtdF(3).rOp(i).rOp().cwise() *= 2.0;
         this->rEtdF(3).rOp(i).rOp() += this->etdF(1).op(i).op()*this->etdF(2).op(i).op() + this->etdF(2).op(i).op();
         this->rEtdF(3).rOp(i).rOp().cwise() *= 0.125;
      }
   }

   template <typename TSim, int TNOps> void ETDNOperators<TSim, TNOps>::computeSsquaredF3()
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
