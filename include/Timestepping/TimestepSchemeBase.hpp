/** \file TimestepSchemeBase.hpp
 *  \brief Base for any implementation of timestep scheme
 */

#ifndef TIMESTEPSCHEMEBASE_HPP
#define TIMESTEPSCHEMEBASE_HPP

// Configuration includes
//
#include "Config/SmartPointer.h"
#include "Config/SimulationConfig.hpp"

// System includes
//

// External includes
//

// Project includes
//
#include "Domain/Truncation.hpp"
#include "Timestepping/TimestepParameters.hpp"

namespace EPMDynamo {

   /**
    * \brief Base for any implementation of timestep scheme
    */
   class TimestepSchemeBase
   {
      public:
         /// Typedef from Simulation trait to local scalar type
         typedef SimulationConfig::NumericalScheme::ScalarType    ScalarType;

         /// Typedef from Simulation trait to local scalar type
         typedef EPMSHARED_PTR<ScalarType>  SmartScalarType;

         /**
          * @brief Constructor
          *
          * @param tsteps Timestep parameters
          * @param pTrunc Truncation information
          * @param hasL0 Is the l=0 mode required?
          */
         TimestepSchemeBase(TimestepParameters &tsteps, SmartTruncation pTrunc, bool hasL0);

         /**
          * @brief Destructor
          */
         virtual ~TimestepSchemeBase() {};
         
      protected:
         /**
          * @brief Get the timestep parameters
          */
         TimestepParameters& rTSParams();

         /**
          * @brief Get pointer to variable from previous timestep
          */
         SmartScalarType pOldVar() const;

         /**
          * @brief Get pointer to nTerms from previous timestep
          */
         SmartScalarType pOldNTerms() const;

         /**
          * @brief Get values of variable from previous timestep
          */
         const ScalarType& oldVar() const;

         /**
          * @brief Get values of nTerms from previous timestep
          */
         const ScalarType& oldNTerms() const;

         /**
          * @brief Store values from previous timestep
          *
          * @param var Old variable values
          * @param nTerms Old non linear terms values
          */
         void storeOld(const ScalarType& var, const ScalarType& nTerms);

         /**
          * @brief Restore values from previous timestep
          *
          * @param rVar Storage to recover old variable value
          * @param nTerms Storage to recover old non linear terms values
          */
         void restoreOld(ScalarType& rVar, ScalarType& nTerms);

      private:
         /**
          * @brief Reference to a Timestep Parameters object
          */
         TimestepParameters&  mrTStepParams;

         /**
          * @brief Values of the variable of previous timestep
          */
         SmartScalarType  mpOldVar;

         /**
          * @brief Values of the nTerms of previous timestep
          */
         SmartScalarType  mpOldNTerms;

         /**
          * @brief Set values of variable from previous timestep
          */
         ScalarType& rOldVar();

         /**
          * @brief Set values of nTerms from previous timestep
          */
         ScalarType& rOldNTerms();
   };

   inline TimestepParameters& TimestepSchemeBase::rTSParams()
   {
      return this->mrTStepParams;
   }

   inline TimestepSchemeBase::SmartScalarType TimestepSchemeBase::pOldVar() const
   {
      return this->mpOldVar;
   }

   inline const TimestepSchemeBase::ScalarType& TimestepSchemeBase::oldVar() const
   {
      return (*this->mpOldVar);
   }

   inline TimestepSchemeBase::ScalarType& TimestepSchemeBase::rOldVar()
   {
      return (*this->mpOldVar);
   }

   inline TimestepSchemeBase::SmartScalarType TimestepSchemeBase::pOldNTerms() const
   {
      return this->mpOldNTerms;
   }

   inline const TimestepSchemeBase::ScalarType& TimestepSchemeBase::oldNTerms() const
   {
      return (*this->mpOldNTerms);
   }

   inline TimestepSchemeBase::ScalarType& TimestepSchemeBase::rOldNTerms()
   {
      return (*this->mpOldNTerms);
   }

}

#endif // TIMESTEPSCHEMEBASE_HPP
