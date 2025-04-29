/** \file XMLReader.cpp 
 *  \brief Source of the implementation of the XML tree reader
 */

// System includes
//

#include <iostream>

// External includes
//

// Class include
//
#include "IO/ASCII/XMLReader.hpp"

// Project includes
//
#include "General/EPMException.hpp"
#include "IO/ASCII/XMLReaderDefs.hpp"

namespace EPMDynamo {

   XMLReader::XMLReader(std::string type, std::string name, std::string ext, std::string header, std::string version)
      : ASCIIReader(name, ext, header, version), mType(type), mContent(), mXML()
   {
   }

   void XMLReader::init()
   {
      if(this->doesIO())
      {
         this->open();

         this->parse();

         this->checkCompatibility();
      }
   }

   void XMLReader::parse()
   {
      if(this->doesIO())
      {
         std::string tmp;
         std::string tmpAll;

         while(mFile.good())
         {
            getline(this->mFile, tmp);
            tmpAll += tmp;
         }

         this->mContent = std::vector<char>(tmpAll.begin(), tmpAll.end());
         this->mContent.push_back('\0');

         this->mXML.parse<0>(&this->mContent[0]);
      }
   }

   void XMLReader::checkCompatibility()
   {
      if(this->doesIO())
      {
         rapidxml::xml_node<> *node = this->mXML.first_node(XMLReaderDefs::GENERALXML.c_str());

         rapidxml::xml_node<> *vnode = node->first_node(XMLReaderDefs::HEADERXML.c_str());
         if(vnode)
         {
            std::string fileHead = vnode->value();
            // Check header
            if(fileHead == this->mHeader)
            {
               vnode = node->first_node(XMLReaderDefs::VERSIONXML.c_str());
               // Read version if present
               if(vnode)
               {
                  std::string fileVers = vnode->value();
                  // Check version
                  if(fileVers == this->mVersion)
                  {
                     vnode = node->first_node(XMLReaderDefs::TYPEXML.c_str());
                     // Read type if present
                     if(vnode)
                     {
                        std::string fileType = vnode->value();
                        // Check type
                        if(fileType == this->mType)
                        {
                        } else
                        {
                           throw EPMException("XMLReader::checkCompatibility", "Wrong file type!");
                        }
                     } else
                     {
                        throw EPMException("XMLReader::checkCompatibility", "Missing file type!");
                     }
                  } else
                  {
                     throw EPMException("XMLReader::checkCompatibility", "Wrong file version!");
                  }
               } else
               {
                  throw EPMException("XMLReader::checkCompatibility", "Missing file version!");
               }
            } else
            {
               throw EPMException("XMLReader::checkCompatibility", "Wrong file header!");
            }
         } else
         {
            throw EPMException("XMLReader::checkCompatibility", "Missing file header!");
         }
      }
   }

}
