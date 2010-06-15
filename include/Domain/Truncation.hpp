/** \file Truncation.hpp
 *  \brief Complete truncation information (global, local and core)
 */

#ifndef TRUNCATION_HPP
#define TRUNCATION_HPP

// Configuration includes
// 
#include "Config/Parallelisation.h"
#include "Config/SmartPointer.h"

// System includes
//

// External includes
//

// Project includes
//
#include "Domain/SimulationTruncation.hpp"
#include "Domain/CoreTruncation.hpp"
#include "Domain/LocalTruncation.hpp"

namespace EPMDynamo {

   /**
    * @brief This clas contains all required truncation information
    */
   class Truncation
   {
      public:
         /**
         * @brief Constructor
         *
         * @param maxN  Highest order index for radial dimension
         * @param nR    Number of radial grid points
         * @param maxL  Highest Harmonic degree
         * @param maxM  Highest Harmonic order
         * @param mp    Longitudinal symmetry assumption
         * @param nCore Number of cores used
         */
         Truncation(const int maxN, const int nR, const int maxL, const int maxM, const int mp, const int nCore);

         /**
         * @brief Simple empty destructor
         */
         virtual ~Truncation() {};

         /**
          * @brief Parallelisation information
          */
         WorkflowType& para();
         
         /**
          * @brief Get the general Simulation wide truncation information
          */
         const SmartSimTrunc  sim() const;
         
         /**
          * @brief Get the local core truncation information
          */
         const SmartLocalTrunc  local() const;
         
         /**
          * @brief Get the remote cores truncation information
          */
         const SmartCoreTrunc  remote(const int id) const;

         /**
          * @brief Set the radial grid (isn't know at creation time)
          *
          * @param rad Radial grid
          */
         void setRadialGrid(SmartArray rad);

         /**
          * @brief Set the theta grid (isn't know at creation time)
          *
          * @param theta Theta grid
          * @param cTh  Cos Theta grid
          * @param sTh Sin Theta grid
          */
         void setThetaGrid(SmartArray theta, SmartArray cTh, SmartArray sTh);
         
      protected:

      private:
         /**
          * @brief Parallelisation information and control
          */
         WorkflowType mPara;

         /**
          * @brief Simulation wide truncation information
          */
         SmartSimTrunc mpSimulation;

         /**
          * @brief Local core truncation information
          */
         SmartLocalTrunc   mpLocal;

         /**
          * @brief A vector of core truncation information
          */
         std::vector<SmartCoreTrunc>   mRemote;

         /**
          * @brief Split the workload
          *
          * \bug Requires correct implemenation of the load splitters
          */
         void splitLoad();
   };

   inline WorkflowType& Truncation::para()
   {
      return this->mPara;
   }

   inline const SmartSimTrunc Truncation::sim() const
   {
      return this->mpSimulation;
   }

   inline const SmartLocalTrunc Truncation::local() const
   {
      return this->mpLocal;
   }

   inline const SmartCoreTrunc Truncation::remote(const int id) const
   {
      return this->mRemote.at(id);
   }

   inline void Truncation::setRadialGrid(SmartArray rad)
   {
      // Set simulation wide radial grid
      this->sim()->rad()->setRadGrid(rad);

      // Set local core radial grid
      this->local()->rtp()->setRadialGrid(rad);
   }

   inline void Truncation::setThetaGrid(SmartArray theta, SmartArray cTh, SmartArray sTh)
   {
      // Set simulation wide theta grid
      this->sim()->hoz()->setThetaGrid(theta, cTh, sTh);

      // Set local core theta grid
      this->local()->rtp()->setThetaGrid(cTh, sTh);
   }

   /// Typedef for a smart Truncation
   typedef EPMSHARED_PTR<Truncation> SmartTruncation;
}

#endif // TRUNCATION_HPP
