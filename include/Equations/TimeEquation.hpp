/** \file TimeEquation.hpp
 *  \brief General implementation of a time dependend equation
 */

#ifndef TIMEEQUATION_HPP
#define TIMEEQUATION_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "General/EPMTypedefs.hpp"
#include "Timestepping/TimestepParameters.hpp"
#include "PhysicalFields/PhysicalVariable.hpp"

namespace EPMDynamo {

   /**
    * @brief General implementation of a time steppable equation.
    *
    * It provides the minimal required methods to use it in a simulation
    *
    * \tparam TSimTypeType Type of the simulation
    * \tparam TFieldType Type of the field (scalar, toroidal/poloidal, QST)
    */
   template <typename TSimType, typename TFieldType > class TimeEquation
   {
      public:
         /// Typedef from Simulation trait to local transform type
         typedef typename SimulationTraits<TSimType>::TransformType    TransformType;

         /**
          *
          * \brief Constructor simply sets the reference to the unknown and the transform
          *
          * \param unknown Reference to the unknown variable of equation
          * \param transform  Reference to the transform object
          * \param tsteps Timestep parameters
          *
          */
         TimeEquation(TFieldType &unknown, TransformType &transform, TimestepParameters& tsteps);

         /**
          * @brief Simple empty destructor
          */
         virtual ~TimeEquation() {};

         /**
          * @brief Pure virtual method for updating the equation RTP values
          *
          * \param step Current step in a multistep transform
          */
         virtual void updateRTP(const int step) = 0;

         /**
          * @brief Pure virtual method for updating the equation RHS values
          */
         virtual void updateRHS() = 0;

         /**
          * @brief Pure virtual method for transforming the equation values to spectral space
          *
          * \param step Current step in a multistep transform
          */
         virtual void transformRHS(const int step) = 0;

         /**
          * @brief Pure virtual method for timestepping the equation
          */
         virtual void timestep() = 0;
         
      protected:

         /**
          * @brief Boundary condition counter
          */
         int mBCCounter;

         /**
          * @brief Reference variable to the unknown of the equation
          */
         TFieldType&   mrX;

         /**
          * @brief Variable for non linear computations
          */
         TFieldType    mNTerms;

         /**
          * @brief Reference to the transform object
          *
          * This is a reference to a TTransform object which should allow to pass
          * in all the currently implemented and considered transforms
          */
         TransformType& mrTransform;

         /**
          * @brief Reference to the timestep parameters
          */
         TimestepParameters&  mrTStepParams;

         /**
          * @brief Check if all BCs are seet
          */
         bool hasAllBCs() const;

      private:
   };

   template <typename TSimType, typename TFieldType> TimeEquation<TSimType, TFieldType>::TimeEquation(TFieldType &unknown, typename TimeEquation<TSimType, TFieldType>::TransformType &transform, TimestepParameters& tsteps)
      : mBCCounter(-1), mrX(unknown), mNTerms(mrX.oc().trunc(), transform), mrTransform(transform), mrTStepParams(tsteps)
   {
   }

   template < typename TSimType, typename TFieldType> bool TimeEquation<TSimType, TFieldType>::hasAllBCs() const
   {
      // Check that the boundary conditions counter reached zero
      if(this->mBCCounter == 0)
      {
         return true;
      }
      else
      {
         return false;
      }
   }
}

#endif // TIMEEQUATION_HPP
