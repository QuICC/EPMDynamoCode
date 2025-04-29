/** \file HDF5NWriter.hpp
 *  \brief Implementation of a numbering HDF5 file
 */

#ifndef HDF5NWRITER_HPP
#define HDF5NWRITER_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "General/EPMTypedefs.hpp"
#include "IO/HDF5/HDF5Writer.hpp"

namespace EPMDynamo {

   /**
    * @brief Class describes an hdf5 file writer where the data gets written to a new file every time
    *          a number is appended to the name.
    */
   class HDF5NWriter: public HDF5Writer
   {
      public:
         /**
         * @brief Constructor
         *
         * @param name File name
         * @param ext File extension
         * @param header File header
         * @param version Version string
         */
         HDF5NWriter(std::string name, std::string ext, std::string header, std::string version);

         /**
         * @brief Destructor
         */
         virtual ~HDF5NWriter() {};

         /**
          * @brief Initialise the file
          */
         void init();

         /**
          * @brief Write the content
          */
         virtual void write() = 0;

         /**
          * @brief Finalise the file
          */
         virtual void finalise();

         /**
          * @brief Change the basename
          *
          * @param base New basename
          */
         void changeBasename(std::string base);
         
      protected:

         /**
          * @brief Update the filename with new number
          */
         void updateName();

         /**
          * @brief Operation to perform just before writing
          */
         void preWrite();

         /**
          * @brief Operations to perform just after writing
          */
         void postWrite();

      private:

         /**
          * @brief Counter for the file number
          */
         int mCounter;

         /**
          * @brief Base name used for appending the number
          */
         std::string mBaseName;

         /**
          * @brief Width of the zero fill for the file id
          */
         static const int msIDWidth;
   };

}

#endif // HDF5WRITER_HPP
