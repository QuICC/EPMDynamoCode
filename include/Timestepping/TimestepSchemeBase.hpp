/** \file TimestepSchemeBase.hpp
 *  \brief Base for any implementation of timestep scheme
 */

#ifndef TIMESTEPSCHEMEBASE_HPP
#define TIMESTEPSCHEMEBASE_HPP

// Configuration includes
//
#include "Config/SmartPointer.h"

// System includes
//

// External includes
//

// Project includes
//
#include "Domain/Truncation.hpp"
#include "Timestepping/TimestepParameters.hpp"
#include "Simulations/Traits/SimulationTraits.hpp"

namespace EPMDynamo {

   /**
    * \brief Base for any implementation of timestep scheme
    *
    * \tparam TSimType Type of the simulation
    */
   template <typename TSimType> class TimestepSchemeBase
   {
      public:
         /// Typedef from Simulation trait to local scalar type
         typedef typename TSimType::ScalarType    ScalarType;

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

   template <typename TSimType> inline TimestepParameters& TimestepSchemeBase<TSimType>::rTSParams()
   {
      return this->mrTStepParams;
   }

   template <typename TSimType> inline typename TimestepSchemeBase<TSimType>::SmartScalarType TimestepSchemeBase<TSimType>::pOldVar() const
   {
      return this->mpOldVar;
   }

   template <typename TSimType> inline const typename TimestepSchemeBase<TSimType>::ScalarType& TimestepSchemeBase<TSimType>::oldVar() const
   {
      return (*this->mpOldVar);
   }

   template <typename TSimType> inline typename TimestepSchemeBase<TSimType>::ScalarType& TimestepSchemeBase<TSimType>::rOldVar()
   {
      return (*this->mpOldVar);
   }

   template <typename TSimType> inline typename TimestepSchemeBase<TSimType>::SmartScalarType TimestepSchemeBase<TSimType>::pOldNTerms() const
   {
      return this->mpOldNTerms;
   }

   template <typename TSimType> inline const typename TimestepSchemeBase<TSimType>::ScalarType& TimestepSchemeBase<TSimType>::oldNTerms() const
   {
      return (*this->mpOldNTerms);
   }

   template <typename TSimType> inline typename TimestepSchemeBase<TSimType>::ScalarType& TimestepSchemeBase<TSimType>::rOldNTerms()
   {
      return (*this->mpOldNTerms);
   }

   template <typename TSimType> TimestepSchemeBase<TSimType>::TimestepSchemeBase(TimestepParameters &tsteps, SmartTruncation pTrunc, bool hasL0)
      : mrTStepParams(tsteps), mpOldVar(new ScalarType(pTrunc, hasL0)), mpOldNTerms(new ScalarType(pTrunc, hasL0)) 
   {
   }

   template <typename TSimType> void TimestepSchemeBase<TSimType>::storeOld(const typename TimestepSchemeBase<TSimType>::ScalarType& var, const typename TimestepSchemeBase<TSimType>::ScalarType& nTerms)
   {
      // Get number of harmonic degrees
      int nL = this->oldNTerms().nL();
      // Get minimum harmonic degrees
      const int l0 = this->oldNTerms().minL();

      // Loop over degrees
      for(int l = l0; l < nL; ++l)
      {
         // Store previous timestep variables
         this->rOldVar().rLShell(l) = var.lshell(l);

         // Store previous timestep nterms
         this->rOldNTerms().rLShell(l) = nTerms.lshell(l);
      }
   }

   template <typename TSimType> void TimestepSchemeBase<TSimType>::restoreOld(typename TimestepSchemeBase<TSimType>::ScalarType& rVar, typename TimestepSchemeBase<TSimType>::ScalarType& nTerms)
   {
      // Get number of harmonic degrees
      int nL = this->oldNTerms().nL();
      // Get minimu harmonic degrees
      const int l0 = this->oldNTerms().minL();

      // Loop over degrees
      for(int l = l0; l < nL; ++l)
      {
         // Recover previous timestep variables
         rVar.rLShell(l) = this->oldVar().lshell(l);

         // Recover previous timestep nterms
         nTerms.rLShell(l) = this->oldNTerms().lshell(l);
      }
   }
}

#endif // TIMESTEPSCHEMEBASE_HPP
