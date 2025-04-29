/** \file XMLReader.hpp
 *  \brief Implementation of an XML reader
 */

#ifndef XMLREADER_HPP
#define XMLREADER_HPP

// System includes
//
#include <fstream>
#include <sstream>

// External includes
//
#include "rapidxml.hpp"

// Project includes
//
#include "General/EPMException.hpp"
#include "IO/ASCII/ASCIIReader.hpp"

namespace EPMDynamo {

   /**
    * @brief Simpler reader class for XML data (based on rapidXML)
    *
    * This will mainly be used for reading in the parameters
    */
   class XMLReader: public ASCIIReader
   {
      public:
         /**
         * @brief Constructor 
         *
         * @param type Type of the file
         * @param name File name
         * @param ext Extension of the file
         * @param header Header of the file
         * @param version Version string of the file
         */
         XMLReader(std::string type, std::string name, std::string ext, std::string header, std::string version);

         /**
         * @brief Destructor
         */
         virtual ~XMLReader() {};

         /**
          * @brief Initialise the file
          */
         void init();

         /**
          * @brief Read the content
          */

         virtual void read() = 0;
         
      protected:
         /**
          * @brief Type of the file
          */
         std::string mType;

         /**
          * @brief Content of file converted into a vector
          */
         std::vector<char> mContent;

         /**
          * @brief XML interface
          */
         rapidxml::xml_document<>   mXML;

         /**
          * @brief Parse the XML content
          */
         void parse();

         /**
          * @brief Check compatibility of opened file
          */
         void checkCompatibility();

         /**
          * @brief templated function to read value from XML tree
          *
          * @param val Storage for read value
          * @param base XML tree node
          * @param tag XML tag to get value
          *
          * \tparam T Type of the value to read
          */
         template <typename T> void readValue(T &val, rapidxml::xml_node<> *base, const std::string &tag);

      private:
   };

   template <typename T>  void XMLReader::readValue(T &val, rapidxml::xml_node<> *base, const std::string &tag)
   {
      std::istringstream   iss;
      std::string tmp;

      rapidxml::xml_node<> *vnode = base->first_node(tag.c_str());
      if(vnode)
      {
         tmp = vnode->value();
         iss.str(tmp);
         iss >> val;
         iss.clear();
      } else
      {
         throw EPMException("XMLReader::readValue", "Unknown tag!");
      }
   }

}

#endif // XMLREADER_HPP
