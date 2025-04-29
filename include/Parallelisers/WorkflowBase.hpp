/** \file WorkflowBase.hpp
 *  \brief Base of the implementation of the workflow control
 */

#ifndef WORKFLOWBASE_HPP
#define WORKFLOWBASE_HPP

// System includes
//

// External includes
//

// Project includes
//

namespace EPMDynamo {

   /**
    * @brief This class give basic information required for setting up the workflow
    */
   class WorkflowBase
   {
      public:
         /**
          * @brief Constructor
          *
          * @param nCore Number of cores
          */
         WorkflowBase(const int nCore);

         /**
          * @brief Destructor
          */
         virtual ~WorkflowBase() {};

         /**
          * @brief Get the number of cores
          */
         int nCore() const;

         /**
          * @brief Get the MPI ID
          */
         int id() const;
         
      protected:
         /**
          * @brief Number of cores
          */
         int mNcore;

         /**
          * @brief ID of the current core (ie MPI rank)
          */
         int mCoreID;

         /**
          * @brief Check compatibilty of number of cores
          *
          * @param cores Number of cores that have been requested
          */
         void checkCores(const int cores);

      private:
   };

   inline int WorkflowBase::nCore() const
   {
      return this->mNcore;
   }

   inline int WorkflowBase::id() const
   {
      return this->mCoreID;
   }

}

#endif // WORKFLOWBASE_HPP
