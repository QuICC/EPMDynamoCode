/** \file ETDSchemeTraits.hpp
 */

#ifndef ETDSCHEMETRAITS_HPP
#define ETDSCHEMETRAITS_HPP

// System includes
//

// External includes
//

// Project includes
//

namespace EPMDynamo {

   // Forward declaration
   template <typename, template <typename> class> class ETDScheme;
   template <typename> class ETD2RKMethod;
   template <typename> class ETD2RKMethod;

   /**
    * @brief ETD traits class
    */
   template <typename TSim> class ETDSchemeTraits
   {
      public:
         /// Typedef for the timestepper without poisson step
         typedef ETDScheme<TSim, ETD2RKMethod>   Timestepper;  

         /// Typedef for the timestepper with poisson step
         typedef ETDScheme<TSim, ETD2RKMethod>   PoissonTimestepper;  

         /// Typedef for the timestepper with poisson step
         typedef ETDTimestepControl   TimestepControl;
         
      protected:

      private:
         /**
          * @brief Constructor
          */
         ETDSchemeTraits();

         /**
          * @brief Destructor
          */
         virtual ~ETDSchemeTraits();
   };

}

#endif // ETDSCHEMETRAITS_HPP
