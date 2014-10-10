/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2013 University of Muenster, Germany.                        *
 * Visualization and Computer Graphics Group <http://viscg.uni-muenster.de>        *
 * For a list of authors please refer to the file "CREDITS.txt".                   *
 *                                                                                 *
 * This file is part of the Voreen software package. Voreen is free software:      *
 * you can redistribute it and/or modify it under the terms of the GNU General     *
 * Public License version 2 as published by the Free Software Foundation.          *
 *                                                                                 *
 * Voreen is distributed in the hope that it will be useful, but WITHOUT ANY       *
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR   *
 * A PARTICULAR PURPOSE. See the GNU General Public License for more details.      *
 *                                                                                 *
 * You should have received a copy of the GNU General Public License in the file   *
 * "LICENSE.txt" along with this file. If not, see <http://www.gnu.org/licenses/>. *
 *                                                                                 *
 * For non-commercial academic use see the license exception specified in the file *
 * "LICENSE-academic.txt". To get information about commercial licensing please    *
 * contact the authors.                                                            *
 *                                                                                 *
 ***********************************************************************************/

#include "voreen/core/io/serialization/xmldeserializer.h"
#include "voreen/core/voreenapplication.h"
#include "voreen/core/voreenmodule.h"
#include "voreen/core/animation/animation.h"

namespace voreen {

const std::string XmlDeserializer::loggerCat_ = "voreen.XmlDeserializer";

XmlDeserializer::XmlDeserializer(std::string documentPath)
    : XmlSerializerBase()
    , documentPath_(documentPath)
{
    // register application (as proxy for modules)
    if (VoreenApplication::app()) {
        registerFactory(VoreenApplication::app());
        registerFactories(VoreenApplication::app()->getSerializerFactories());
    }
    else {
        LWARNING("VoreenApplication not instantiated");
    }
}

XmlDeserializer::~XmlDeserializer() {
    for (UnresolvedReferenceMapType::iterator mapIt = unresolvedReferenceMap_.begin();
        mapIt != unresolvedReferenceMap_.end(); ++mapIt)
    {
        for (ReferenceResolverListType::iterator listIt = mapIt->second.begin();
            listIt != mapIt->second.end(); ++listIt)
        {
            delete *listIt;
        }
    }
}

std::string XmlDeserializer::getDocumentPath() const {
    return documentPath_;
}

void XmlDeserializer::deserialize(const std::string& key, bool& data)
    throw (SerializationException)
{
    std::string boolValue;
    deserializeSimpleTypes(key, boolValue);

    std::transform(boolValue.begin(), boolValue.end(), boolValue.begin(), tolower);

    if (boolValue == "true" || boolValue == "1")
        data = true;
    else if (boolValue == "false" || boolValue == "0")
        data = false;
    else
        raise(XmlSerializationFormatException("XML node with key '" + key + "' contains unknown bool value."));
}

void XmlDeserializer::deserialize(const std::string& key, char& data)
    throw (SerializationException)
{
    deserializeSimpleTypes(key, data);
}

void XmlDeserializer::deserialize(const std::string& key, signed char& data)
    throw (SerializationException)
{
    deserializeSimpleTypes(key, data);
}

void XmlDeserializer::deserialize(const std::string& key, unsigned char& data)
    throw (SerializationException)
{
    deserializeSimpleTypes(key, data);
}

void XmlDeserializer::deserialize(const std::string& key, uint16_t& data)
    throw (SerializationException)
{
    deserializeSimpleTypes(key, data);
}

void XmlDeserializer::deserialize(const std::string& key, int16_t& data)
    throw (SerializationException)
{
    deserializeSimpleTypes(key, data);
}

void XmlDeserializer::deserialize(const std::string& key, uint32_t& data)
    throw (SerializationException)
{
    deserializeSimpleTypes(key, data);
}

#ifdef __APPLE__
void XmlDeserializer::deserialize(const std::string& key, long unsigned int& data)
    throw (SerializationException)
{
    deserializeSimpleTypes(key, data);
}
#endif

void XmlDeserializer::deserialize(const std::string& key, int32_t& data)
    throw (SerializationException)
{
    deserializeSimpleTypes(key, data);
}

void XmlDeserializer::deserialize(const std::string& key, uint64_t& data)
    throw (SerializationException)
{
    deserializeSimpleTypes(key, data);
}

void XmlDeserializer::deserialize(const std::string& key, int64_t& data)
    throw (SerializationException)
{
    deserializeSimpleTypes(key, data);
}

void XmlDeserializer::deserialize(const std::string& key, float& data)
    throw (SerializationException)
{
    deserializeSimpleTypes(key, data);
}

void XmlDeserializer::deserialize(const std::string& key, double& data)
    throw (SerializationException)
{
    deserializeSimpleTypes(key, data);
}

void XmlDeserializer::deserialize(const std::string& key, long double& data)
    throw (SerializationException)
{
    deserializeSimpleTypes(key, data);
}

void XmlDeserializer::deserialize(const std::string& key, std::string& data)
    throw (SerializationException)
{
    deserializeSimpleTypes(key, data);
}

void XmlDeserializer::deserialize(const std::string& key, tgt::vec2& data)
    throw (SerializationException)
{
    deserializeTgtVector(key, data);
}

void XmlDeserializer::deserialize(const std::string& key, tgt::vec3& data)
    throw (SerializationException)
{
    deserializeTgtVector(key, data);
}

void XmlDeserializer::deserialize(const std::string& key, tgt::vec4& data)
    throw (SerializationException)
{
    deserializeTgtVector(key, data);
}

void XmlDeserializer::deserialize(const std::string& key, tgt::dvec2& data)
    throw (SerializationException)
{
    deserializeTgtVector(key, data);
}

void XmlDeserializer::deserialize(const std::string& key, tgt::dvec3& data)
    throw (SerializationException)
{
    deserializeTgtVector(key, data);
}

void XmlDeserializer::deserialize(const std::string& key, tgt::dvec4& data)
    throw (SerializationException)
{
    deserializeTgtVector(key, data);
}


void XmlDeserializer::deserialize(const std::string& key, tgt::ivec2& data)
    throw (SerializationException)
{
    deserializeTgtVector(key, data);
}

void XmlDeserializer::deserialize(const std::string& key, tgt::ivec3& data)
    throw (SerializationException)
{
    deserializeTgtVector(key, data);
}

void XmlDeserializer::deserialize(const std::string& key, tgt::ivec4& data)
    throw (SerializationException)
{
    deserializeTgtVector(key, data);
}

void XmlDeserializer::deserialize(const std::string& key, tgt::col3& data)
    throw (SerializationException)
{
    deserializeTgtVector(key, data, true);
}

void XmlDeserializer::deserialize(const std::string& key, tgt::col4& data)
    throw (SerializationException)
{
    deserializeTgtVector(key, data, true);
}

void XmlDeserializer::deserialize(const std::string& key, tgt::mat2& data)
    throw (SerializationException)
{
    tgt::vec2 row0, row1;
    deserializeTgtVector(key+".row0", row0);
    deserializeTgtVector(key+".row1", row1);
    data = tgt::mat2(row0, row1);
}

void XmlDeserializer::deserialize(const std::string& key, tgt::mat3& data)
    throw (SerializationException)
{
    tgt::vec3 row0, row1, row2;
    deserializeTgtVector(key+".row0", row0);
    deserializeTgtVector(key+".row1", row1);
    deserializeTgtVector(key+".row2", row2);
    data = tgt::mat3(row0, row1, row2);
}

void XmlDeserializer::deserialize(const std::string& key, tgt::mat4& data)
    throw (SerializationException)
{
    tgt::vec4 row0, row1, row2, row3;
    deserializeTgtVector(key+".row0", row0);
    deserializeTgtVector(key+".row1", row1);
    deserializeTgtVector(key+".row2", row2);
    deserializeTgtVector(key+".row3", row3);
    data = tgt::mat4(row0, row1, row2, row3);
}

void XmlDeserializer::deserialize(const std::string& key, tgt::Matrix2d& data)
    throw (SerializationException)
{
    tgt::dvec2 row0, row1;
    deserializeTgtVector(key+".row0", row0);
    deserializeTgtVector(key+".row1", row1);
    data = tgt::Matrix2d(row0, row1);
}

void XmlDeserializer::deserialize(const std::string& key, tgt::Matrix3d& data)
    throw (SerializationException)
{
    tgt::dvec3 row0, row1, row2;
    deserializeTgtVector(key+".row0", row0);
    deserializeTgtVector(key+".row1", row1);
    deserializeTgtVector(key+".row2", row2);
    data = tgt::Matrix3d(row0, row1, row2);
}

void XmlDeserializer::deserialize(const std::string& key, tgt::Matrix4d& data)
    throw (SerializationException)
{
    tgt::dvec4 row0, row1, row2, row3;
    deserializeTgtVector(key+".row0", row0);
    deserializeTgtVector(key+".row1", row1);
    deserializeTgtVector(key+".row2", row2);
    deserializeTgtVector(key+".row3", row3);
    data = tgt::Matrix4d(row0, row1, row2, row3);
}

void XmlDeserializer::deserialize(const std::string& key, Serializable& data)
    throw (SerializationException)
{
    deserializeSimpleTypes(key, data);
}

TiXmlElement* XmlDeserializer::getNextXmlElement(const std::string& key)
    throw (SerializationException)
{
    TiXmlElement* element = node_->FirstChildElement(key);
    while (element) {
        // Was node not visited before?
        if (visitedNodes_.find(element) == visitedNodes_.end())
        {
            visitedNodes_.insert(element);
            return element;
        }

        element = element->NextSiblingElement(key);
    }

    raise(XmlSerializationNoSuchDataException("No further XML node with key '" + key + "' found."));
    return 0;
}

void XmlDeserializer::findUnresolvableReferences(
    TiXmlElement* node,
    ReferenceIdSetType& references,
    ReferenceIdSetType& resolvableReferences) const
{
    const std::string* id = node->Attribute(XmlSerializationConstants::IDATTRIBUTE);
    const std::string* ref = node->Attribute(XmlSerializationConstants::REFERENCEATTRIBUTE);

    if (id)
        resolvableReferences.insert(*id);

    if (ref)
        references.insert(*ref);

    for (TiXmlElement* child = node->FirstChildElement(); child != 0; child = child->NextSiblingElement())
        findUnresolvableReferences(child, references, resolvableReferences);
}

XmlDeserializer::ReferenceIdListType XmlDeserializer::findUnresolvableReferences() const
{
    ReferenceIdListType unresolvableReferences;

    ReferenceIdSetType references;
    ReferenceIdSetType resolvableReferences;
    findUnresolvableReferences(node_->ToElement(), references, resolvableReferences);

    for (ReferenceIdSetType::const_iterator it = references.begin(); it != references.end(); ++it)
        if (resolvableReferences.find(*it) == resolvableReferences.end())
            unresolvableReferences.push_back(*it);

    return unresolvableReferences;
}

void XmlDeserializer::read(std::istream& stream, XmlProcessor* xmlProcessor)
    throw (SerializationException)
{
    // Read input stream...
    std::stringbuf buffer;
    do
    {
        // Use 0 character instead of '\n' to minimize the number of get-calls...
        stream.get(buffer, 0);
    } while (stream.good() && !stream.eof()
        && (buffer.sputc(stream.get()) != std::stringbuf::traits_type::eof()));

    // Parse input...
    document_.Parse(buffer.str().c_str());

    TiXmlElement* root = document_.RootElement();

    // Is there no root element?
    if (!root)
        raise(XmlSerializationFormatException(std::string("No root node found.")));

    // Has root node incorrect name?
    if (root->ValueStr() != XmlSerializationConstants::ROOTNODE) {
        raise(XmlSerializationFormatException("XML root node name is '" + root->ValueStr()
            + "' instead of '" + XmlSerializationConstants::ROOTNODE + "'."));
    }

    const std::string* version = root->Attribute(XmlSerializationConstants::VERSIONATTRIBUTE);
    // Is serialization version not set?
    if (!version)
        raise(XmlSerializationFormatException("XML root node has no version attribute."));
    // Does XmlSerializer and XmlDeserializer version not match the XML document version?
    if (*version != XmlSerializationConstants::VERSION) {
        raise(XmlSerializationVersionMismatchException("XML document has version " + *version
            + " instead of " + XmlSerializationConstants::VERSION + "."));
    }

    node_ = root;

    // Apply preprocessor, if a XML preprocessor is given...
    if (xmlProcessor)
        xmlProcessor->process(document_);

    std::vector<std::string> unresolvableReferences = findUnresolvableReferences();
    if (!unresolvableReferences.empty()) {
        std::stringstream idStream;
        for (ReferenceIdListType::iterator it = unresolvableReferences.begin(); it != unresolvableReferences.end(); ++it)
            idStream << (idStream.str().empty() ? "" : ", ") << "'" << *it << "'";

        raise(XmlSerializationReferenceResolvingException(
            "XML document contains the following unresolvable references: " + idStream.str() + "."));
    }
}

void XmlDeserializer::freePointer(void* pointer) {
    for (IdAddressMapType::iterator it = idAddressMap_.begin(); it != idAddressMap_.end(); ++it)
        if (it->second == pointer)
            it->second = 0;
}

void XmlDeserializer::deserializeBinaryBlob(const std::string& key, unsigned char* inputBuffer, size_t reservedMemory)
    throw (SerializationException)
{
    std::string tmp;
    deserialize(key, tmp);
    std::vector<unsigned char> dataVec = base64Decode(tmp);

    if(dataVec.size() != reservedMemory)
        raise(XmlSerializationFormatException("XML node with key '" + key + "' contains a binary blob for which unsufficient memory has been reserved."));
    else {
        //std::copy(inputBuffer, inputBuffer + dataVec.size(), dataVec.begin()); //doesn't work
        memcpy(inputBuffer, &dataVec[0], reservedMemory);
    }
}

void XmlDeserializer::deserializeBinaryBlob(const std::string& key, std::vector<unsigned char>& buffer)
    throw (SerializationException)
{
    std::string tmp;
    deserialize(key, tmp);
    buffer = base64Decode(tmp);
}

std::vector<unsigned char> XmlDeserializer::base64Decode(const std::string& input) {
    //input.erase(std::remove(input.begin(), input.end(), '\n'), input.end());

    const char padCharacter = '=';

    if (input.length() % 4) //Sanity check
            raise(XmlSerializationFormatException("Non-Valid base64!"));

    size_t padding = 0;
    if (input.length()) {
        if (input[input.length()-1] == padCharacter)
                padding++;
        if (input[input.length()-2] == padCharacter)
                padding++;
    }

    //Setup a vector to hold the result
    std::vector<unsigned char> decodedBytes;
    decodedBytes.reserve(((input.length()/4)*3) - padding);
    long temp = 0; //Holds decoded quanta
    std::string::const_iterator cursor = input.begin();
    while (cursor < input.end()) {
        for (size_t quantumPosition = 0; quantumPosition < 4; quantumPosition++) {
            temp <<= 6;
            if       (*cursor >= 0x41 && *cursor <= 0x5A) // This area will need tweaking if
                temp |= *cursor - 0x41;                       // you are using an alternate alphabet
            else if  (*cursor >= 0x61 && *cursor <= 0x7A)
                temp |= *cursor - 0x47;
            else if  (*cursor >= 0x30 && *cursor <= 0x39)
                temp |= *cursor + 0x04;
            else if  (*cursor == 0x2B)
                temp |= 0x3E; //change to 0x2D for URL alphabet
            else if  (*cursor == 0x2F)
                temp |= 0x3F; //change to 0x5F for URL alphabet
            else if  (*cursor == padCharacter) { //pad
                switch( input.end() - cursor ) {
                    case 1: //One pad character
                        decodedBytes.push_back((temp >> 16) & 0x000000FF);
                        decodedBytes.push_back((temp >> 8 ) & 0x000000FF);
                        return decodedBytes;
                    case 2: //Two pad characters
                        decodedBytes.push_back((temp >> 10) & 0x000000FF);
                        return decodedBytes;
                    default:
                        raise(XmlSerializationFormatException("Invalid Padding in Base 64!"));
                }
            }  else
                raise(XmlSerializationFormatException("Non-Valid Character in Base 64!"));
            cursor++;
        }
        decodedBytes.push_back((temp >> 16) & 0x000000FF);
        decodedBytes.push_back((temp >> 8 ) & 0x000000FF);
        decodedBytes.push_back((temp      ) & 0x000000FF);
    }

    return decodedBytes;
}


} // namespace
