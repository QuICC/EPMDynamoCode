/** \file ASCIIFieldWriter.hpp
 *  \brief Implementation of an ASCII field values writer
 */

#ifndef ASCIIFIELDWRITER_HPP
#define ASCIIFIELDWRITER_HPP

// System includes
//
#include <fstream>

// External includes
//

// Project includes
//

namespace EPMDynamo {

   /**
    * @brief Small ASCII file writing class for field values
    *
    * \tparam TSimType Type of the simulation
    * \tparam TFieldType Type of the field
    * \tparam TWriterType Type of the writer
    */
   template <typename TSimType, typename TFieldType, typename TWriterType> class ASCIIFieldWriter: public TWriterType
   {
      public:
         /**
         * @brief Constructor
         *
         * @param var Field variable
         * @param name Filename
         * @param ext File extension
         * @param header Header string of file
         * @param version  Version string of file 
         */
         ASCIIFieldWriter(TFieldType &var, std::string name, std::string ext, std::string header, std::string version);

         /**
         * @brief Destructor
         */
         virtual ~ASCIIFieldWriter() {};
         
      protected:

         /**
          * @brief Reference to the variable to get data from
          */
         TFieldType& mrVar;

      private:
   };

   template <typename TSimType, typename TFieldType, typename TWriterType> ASCIIFieldWriter<TSimType, TVariable, TWriterType>::ASCIIFieldWriter(TVariable &var, std::string name, std::string ext, std::string header, std::string version)
      : TWriterType(name, ext, header, version), mrVar(var)
   {
   }
   
}

#endif // ASCIIFIELDWRITER_HPP
