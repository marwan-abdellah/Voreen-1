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

#include "voreen/core/datastructures/transfunc/transfunc1dkeys.h"
#include "voreen/core/datastructures/transfunc/transfuncmappingkey.h"

#include "voreen/core/datastructures/transfunc/preintegrationtable.h"

#ifdef VRN_MODULE_DEVIL
    #include <IL/il.h>
#endif

#include "tgt/logmanager.h"

#include <tinyxml/tinyxml.h>
#include <fstream>
#include <sstream>

using tgt::col4;
using tgt::ivec3;
using tgt::ivec4;

namespace voreen {

PreIntegrationTableMap::PreIntegrationTableMap(TransFunc1DKeys* tf, size_t maxSize) : tf_(tf) {
    //maxSize_ must be >= 1
    if (maxSize == 0)
        maxSize_ = 1;
    else
        maxSize_ = maxSize;
}

PreIntegrationTableMap::~PreIntegrationTableMap() {
    clear();
}

void PreIntegrationTableMap::setTransFunc(TransFunc1DKeys* tf) {
    tf_ = tf;
}

void PreIntegrationTableMap::clear() {
    //delete all pre-integration tables, then clear map and queue
    std::map<PreIntegrationTableAttributes, PreIntegrationTable*>::iterator it = tableMap_.begin();
    for (; it!=tableMap_.end(); ++it)
        delete it->second;

    tableMap_.clear();
    queue_ = std::queue<PreIntegrationTableAttributes>();
}

const PreIntegrationTable* PreIntegrationTableMap::getPreIntegrationTable(float samplingStepSize, size_t dimension, bool useIntegral, bool computeOnGPU, tgt::Shader* program) {

    if (!tf_)
        return 0;

    size_t dim;

    if (dimension == 0)
        dim = static_cast<size_t>(std::min(tf_->getDimensions().x, 1024));
    else
        dim = dimension;

    //create key to pre-integration table map
    PreIntegrationTableAttributes attributes;
    attributes.samplingStepSize_ = samplingStepSize;
    attributes.dimension_ = dim;
    attributes.useIntegral_ = useIntegral;
    attributes.computeOnGPU_ = computeOnGPU;

    //try to find the right pre-integration table
    std::map<PreIntegrationTableAttributes, PreIntegrationTable*>::iterator it = tableMap_.find(attributes);
    if (it != tableMap_.end())
        return (it->second);
    else {
        //if size == max: delete the first
        if (tableMap_.size() == maxSize_) {
            it = tableMap_.find(queue_.front());
            delete it->second;
            tableMap_.erase(it);
            queue_.pop();
        }
        //create new pre-integration table and add it to the map
        PreIntegrationTable* table = new PreIntegrationTable(tf_, dim, samplingStepSize, useIntegral, computeOnGPU, program);
        tableMap_.insert(std::make_pair(attributes,table));
        queue_.push(attributes);
        it = tableMap_.find(attributes);
        return (it->second);
    }
}

//-------------------------------------------------------------------------------------------------

// some stream helper functions
inline int readInt(std::ifstream& stream);
inline short readShort(std::ifstream& stream);
inline double readDouble(std::ifstream& stream);

const std::string TransFunc1DKeys::loggerCat_("voreen.TransFunc1DKeys");

TransFunc1DKeys::TransFunc1DKeys(int width)
    : TransFunc(width, 1, 1, GL_RGBA, GL_UNSIGNED_BYTE, tgt::Texture::NEAREST)
    , threshold_(0.f, 1.f)
    , domain_(0.f, 1.f)
    , preIntegrationTableMap_(NULL, 10)
    , preIntegrationProgram_(0)
{
    loadFileFormats_.push_back("tfi");
    loadFileFormats_.push_back("lut");
    loadFileFormats_.push_back("table");
    loadFileFormats_.push_back("plist");
#ifdef VRN_MODULE_DEVIL
    loadFileFormats_.push_back("jpg");
    loadFileFormats_.push_back("png");
#endif

    saveFileFormats_.push_back("tfi");
    saveFileFormats_.push_back("lut");
#ifdef VRN_MODULE_DEVIL
    saveFileFormats_.push_back("png");
#endif

    setToStandardFunc();

    preIntegrationTableMap_.setTransFunc(this);
}

TransFunc1DKeys::TransFunc1DKeys(const TransFunc1DKeys& tf)
    : TransFunc(tf.dimensions_.x, tf.dimensions_.y, tf.dimensions_.z,
                tf.format_, tf.dataType_, tf.filter_)
      , preIntegrationTableMap_(NULL,10)
      , preIntegrationProgram_(0)
{
    updateFrom(tf);
    //pre-integration table map is not copied
    preIntegrationTableMap_.setTransFunc(this);
}

TransFunc1DKeys::~TransFunc1DKeys() {
    for (size_t i = 0; i < keys_.size(); ++i)
        delete keys_[i];

    preIntegrationTableMap_.clear();

    if (preIntegrationProgram_)
        ShdrMgr.dispose(preIntegrationProgram_);
}

bool TransFunc1DKeys::operator==(const TransFunc1DKeys& tf) {
    if(threshold_ != tf.threshold_)
        return false;
    if(domain_ != tf.getDomain(0))
        return false;
    if(gammaValue_ != tf.gammaValue_)
        return false;
    if(alphaMode_ != tf.alphaMode_)
        return false;

    if (keys_.size() != tf.keys_.size())
        return false;
    for (size_t i = 0; i < keys_.size(); ++i) {
        if (*(keys_[i]) != *(tf.keys_[i]))
            return false;
    }

    return true;
}

bool TransFunc1DKeys::operator!=(const TransFunc1DKeys& tf) {
    return !(*this == tf);
}

TransFunc* TransFunc1DKeys::clone() const {
    TransFunc1DKeys* func = new TransFunc1DKeys();

    func->updateFrom(this);

    return func;
}


void TransFunc1DKeys::updateFrom(const TransFunc* transfunc) {
    tgtAssert(transfunc, "null pointer passed");

    const TransFunc1DKeys* tf1DKeys = dynamic_cast<const TransFunc1DKeys*>(transfunc);
    if (!tf1DKeys) {
        LWARNING("updateFrom(): passed parameter is not of type TransFunc1DKeys");
        return;
    }

    dimensions_ = tf1DKeys->getDimensions();
    format_ = tf1DKeys->format_;
    dataType_ = tf1DKeys->dataType_;
    filter_ = tf1DKeys->filter_;
    gammaValue_ = tf1DKeys->gammaValue_;
    alphaMode_ = tf1DKeys->alphaMode_;

    threshold_ = tf1DKeys->threshold_;

    domain_ = tf1DKeys->domain_;

    keys_.clear();
    std::vector<TransFuncMappingKey*>::const_iterator it;
    for (it = tf1DKeys->keys_.begin(); it!=tf1DKeys->keys_.end(); it++) {
        keys_.push_back((*it)->clone());
    }

    textureInvalid_ = true;
}

bool TransFunc1DKeys::isStandardFunc() const {
    if(getDomain() == tgt::vec2(0.0f, 1.0f) && (getNumKeys() == 2)) {
        const TransFuncMappingKey* k0 = getKey(0);
        if((k0->getIntensity() == 0.0f) && !(k0->isSplit()) && (k0->getColorL() == col4(0, 0, 0, 0))) {
            const TransFuncMappingKey* k1 = getKey(1);
            if((k1->getIntensity() == 1.0f) && !(k1->isSplit()) && (k1->getColorL() == col4(255)) &&
                gammaValue_ == 1.f && alphaMode_ == TF_USE_ALPHA)
                return true;
        }
    }
    return false;
}

void TransFunc1DKeys::makeRamp() {
    if(keys_.size() >= 2) {
        float min = keys_[0]->getIntensity();
        float max = keys_[keys_.size()-1]->getIntensity();
        float width = max - min;

        for(size_t i=0; i<keys_.size(); i++) {
            float v = keys_[i]->getIntensity();
            v -= min;
            v /= width;
            keys_[i]->setIntensity(v);
            keys_[i]->setAlphaL(v);
            keys_[i]->setAlphaR(v);
        }
    }
    invalidateTexture();
}

void TransFunc1DKeys::setToStandardFunc() {
    clearKeys();
    keys_.push_back(new TransFuncMappingKey(0.f, col4(0, 0, 0, 0)));
    keys_.push_back(new TransFuncMappingKey(1.f, col4(255)));
    setThresholds(0.f,1.f);
    setGammaValue(1.f);
    setAlphaMode(TF_USE_ALPHA);

    invalidateTexture();
}

tgt::vec4 TransFunc1DKeys::getMeanValue(float segStart, float segEnd) const {
    ivec4 result(0);
    if (!tex_ || textureInvalid_) {
        LWARNING("getMeanValue(): texture is invalid");
        return result;
    }

    float width = static_cast<float>(tex_->getWidth());
    for (int i = static_cast<int>(segStart*width); i < segEnd*width; ++i)
        result += ivec4(tex_->texel<col4>(i));

    return static_cast<tgt::vec4>(result)/(segEnd*width-segStart*width);
}

col4 TransFunc1DKeys::getMappingForValue(float value) const {
    // If there are no keys, any further calculation is meaningless
    if (keys_.empty())
        return col4(0, 0, 0, 0);

    // Restrict value to [0,1]
    value = (value < 0.f) ? 0.f : value;
    value = (value > 1.f) ? 1.f : value;

    //use gamma correction
    if(gammaValue_ != 1.f) {
        value = powf(value,gammaValue_);
    }

    // iterate through all keys until we get to the correct position
    std::vector<TransFuncMappingKey*>::const_iterator keyIterator = keys_.begin();

    while ((keyIterator != keys_.end()) && (value > (*keyIterator)->getIntensity()))
        keyIterator++;

    col4 color;
    if (keyIterator == keys_.begin())
        color = keys_[0]->getColorL();
    else if (keyIterator == keys_.end())
        color = (*(keyIterator-1))->getColorR();
    else{
        // calculate the value weighted by the destination to the next left and right key
        TransFuncMappingKey* leftKey = *(keyIterator-1);
        TransFuncMappingKey* rightKey = *keyIterator;
        float fraction = (value - leftKey->getIntensity()) / (rightKey->getIntensity() - leftKey->getIntensity());
        col4 leftDest = leftKey->getColorR();
        col4 rightDest = rightKey->getColorL();
        color = leftDest;
        color.r += static_cast<uint8_t>((rightDest.r - leftDest.r) * fraction);
        color.g += static_cast<uint8_t>((rightDest.g - leftDest.g) * fraction);
        color.b += static_cast<uint8_t>((rightDest.b - leftDest.b) * fraction);
        color.a += static_cast<uint8_t>((rightDest.a - leftDest.a) * fraction);
    }
    if(alphaMode_ == TF_ONE_ALPHA)
        color.a = 255;
    else if (alphaMode_ == TF_ZERO_ALPHA)
        color.a = 0;

    return color;
}

void TransFunc1DKeys::updateTexture() {

    if (!tex_ || (tex_->getDimensions() != dimensions_))
        createTex();
    tgtAssert(tex_, "No texture");

    //get tresholds of the transfer function
    int front_end = tgt::iround(threshold_.x*dimensions_.x);
    int back_start = tgt::iround(threshold_.y*dimensions_.x);
    //all values before front_end and after back_start are set to zero
    //all other values are taken from the TF
    for (int x = 0; x < front_end; ++x)
        tex_->texel<col4>(x) = col4(0, 0, 0, 0);

    for (int x = front_end; x < back_start; ++x) {
        tex_->texel<col4>(x) = getMappingForValue(static_cast<float>(x) / dimensions_.x);
    }

    for (int x = back_start; x < dimensions_.x; ++x)
        tex_->texel<col4>(x) = col4(0, 0, 0, 0);

    tex_->uploadTexture();
    LGL_ERROR;

    textureInvalid_ = false;
}

void TransFunc1DKeys::setThresholds(const tgt::vec2& thresholds, size_t dimension) {
    threshold_ = thresholds;
    invalidateTexture();
}

tgt::vec2 TransFunc1DKeys::getThresholds(size_t dimension) const {
    return threshold_;
}

int TransFunc1DKeys::getNumKeys() const {
    return static_cast<int>(keys_.size());
}

const TransFuncMappingKey* TransFunc1DKeys::getKey(int i) const {
    return keys_.at(i);
}

TransFuncMappingKey* TransFunc1DKeys::getKey(int i) {
    return keys_.at(i);
}

const std::vector<TransFuncMappingKey*> TransFunc1DKeys::getKeys() const {
    return keys_;
}

void TransFunc1DKeys::setKeys(std::vector<TransFuncMappingKey*> keys) {
    keys_ = keys;
    invalidateTexture();
}

void TransFunc1DKeys::reset(){
    setToStandardFunc();
}

void TransFunc1DKeys::addKey(TransFuncMappingKey* key) {
    // insert key at appropriate location
    if (keys_.empty()) {
        keys_.push_back(key);
        return;
    }
    std::vector<TransFuncMappingKey *>::iterator keyIterator = keys_.begin();
    // Fast-forward to the correct position
    while ((keyIterator != keys_.end()) && (key->getIntensity() > (*keyIterator)->getIntensity()))
        keyIterator++;
    keys_.insert(keyIterator, key);

    invalidateTexture();
}

bool sortFunction(TransFuncMappingKey* a, TransFuncMappingKey* b) {
    return a->getIntensity() < b->getIntensity();
}

void TransFunc1DKeys::updateKey(TransFuncMappingKey* /*key*/) {
    std::sort(keys_.begin(), keys_.end(), sortFunction);

    invalidateTexture();
}

void TransFunc1DKeys::removeKey(TransFuncMappingKey* key) {
    std::vector<TransFuncMappingKey *>::iterator keyIterator = find(keys_.begin(), keys_.end(), key);
    if (keyIterator != keys_.end())
        keys_.erase(keyIterator);
    delete key;

    invalidateTexture();
}

void TransFunc1DKeys::clearKeys() {
    std::vector<TransFuncMappingKey *>::iterator keyIterator = keys_.begin();
    // First delete all the referenced objects in the heap
    while (keyIterator != keys_.end()) {
        delete (*keyIterator);
        ++keyIterator;
    }
    // then delete the entries in the vector
    keys_.clear();

    invalidateTexture();
}

bool TransFunc1DKeys::isEmpty() const {
    return keys_.empty();
}

void TransFunc1DKeys::invalidateTexture() {
    TransFunc::invalidateTexture();
    preIntegrationTableMap_.clear();
}

bool TransFunc1DKeys::save(const std::string& filename) const {
    //look for fileExtension
    std::string fileExtension;
    size_t dotPosition = filename.rfind(".");
    if (dotPosition == std::string::npos)
        return false;
    else
        fileExtension = filename.substr(dotPosition+1);

    if (fileExtension == "tfi")
        return saveTfi(filename);
    else if (fileExtension == "lut")
        return saveLUT(filename);
    else
        return saveImage(filename);
}

bool TransFunc1DKeys::saveTfi(const std::string& filename) const {

    // open file stream
    std::ofstream stream(filename.c_str(), std::ios_base::out);
    if (stream.fail()) {
        LWARNING("Unable to open file " << filename << " for writing.");
        return false;
    }

    // serialize to stream
    bool success = true;
    try {
        XmlSerializer s(filename);
        s.serialize("TransFuncIntensity", this);

        s.write(stream);
        if (stream.bad()) {
            LWARNING("Unable to write to file: " << filename);
            success = false;
        }
        stream.close();
    }
    catch (SerializationException &e) {
        LWARNING("SerializationException: " << e.what());
        stream.close();
        success = false;
    }

    // log result
    if (success)
        LINFO("Saved transfer function to file: " << filename);
    else
        LWARNING("Saving transfer function failed.");

    return success;
}

bool TransFunc1DKeys::saveImage(const std::string& filename) const {
#ifdef VRN_MODULE_DEVIL
    if (!tex_ || textureInvalid_) {
        LWARNING("saveImage(): texture is invalid");
        return false;
    }

    //extract file extension
    std::string fileExtension;
    size_t dotPosition = filename.rfind(".");
    fileExtension = filename.substr(dotPosition+1);

    //IL does _NOT_ overwrite files by default
    ilEnable(IL_FILE_OVERWRITE);
    ILuint img;
    ilGenImages(1, &img);
    ilBindImage(img);

    int width = tex_->getWidth();
    GLubyte* im = new GLubyte[width*4];

    for (int x = 0; x < width; ++x) {
        col4 c = getMappingForValue( static_cast<float>(x)/width );
        im[(x*4)+0] = c.r;
        im[(x*4)+1] = c.g;
        im[(x*4)+2] = c.b;
        im[(x*4)+3] = c.a;
    }

    ilTexImage(width, 1, 1, 4, IL_RGBA, IL_UNSIGNED_BYTE, im);

    if (fileExtension == "png")
        ilSave(IL_PNG, (ILstring)filename.c_str());
    else {
        ilDeleteImages(1, &img);
        delete[] im;
        return false;
    }

    ilDeleteImages(1, &img);
    delete[] im;
    return true;
#else
    LERROR("Saving of " << filename  << " failed: No DevIL support.");
    return false;
#endif // VRN_MODULE_DEVIL
}

bool TransFunc1DKeys::saveLUT(const std::string& filename) const {
    if (!tex_ || textureInvalid_) {
        LWARNING("saveLUT(): texture is invalid");
        return false;
    }

    std::fstream f;
    f.open(filename.c_str(), std::ios::out);
    for (int i = 0; i < tex_->getWidth(); ++i) {
        col4 c = getMappingForValue(static_cast<float>(i)/tex_->getWidth());
        f << static_cast<int>(c.r) << " "
          << static_cast<int>(c.g) << " "
          << static_cast<int>(c.b) << std::endl;
    }

    f.close();
    return true;
}

bool TransFunc1DKeys::load(const std::string& filename) {
    // Extract the file extension
    std::string fileExtension;
    size_t dotPosition = filename.rfind(".");
    if (dotPosition != std::string::npos)
        // => the last (seperating) dot was found
        fileExtension = filename.substr(dotPosition+1);
    else
        return false;

    // switch-statement is not possible, because the compiler doesn't allow a switch by std::string
    // If further extensions shall be added, the pattern should be obvious.
    if (fileExtension == "tfi")
        return loadTfi(filename);
    if (fileExtension == "lut")
        return loadImageJ(filename);
    if (fileExtension == "table")
        return loadTextTable(filename);
    if (fileExtension == "plist")
        return loadOsirixCLUT(filename);
    // "else"
    return loadImage(filename);
}

#ifdef VRN_MODULE_DEVIL
bool TransFunc1DKeys::loadImage(const std::string& filename) {
    ILuint imageID;
    ILboolean success;

    // Load and bind the image
    ilGenImages(1 , &imageID);
    ilBindImage(imageID);
    success = ilLoadImage(const_cast<char*>(filename.c_str()));

    // If the load was unsuccessful => jump out of the method
    if (!success)
        return false;

    ILint imageFormat = ilGetInteger(IL_IMAGE_FORMAT);

    // Convert every color component to unsigned byte value
    success = ilConvertImage(imageFormat , IL_UNSIGNED_BYTE);

    // If the conversion was unsuccessful => jump out of the method
    if (!success)
        return false;

    // Get the width of the image file
    dimensions_.x = ilGetInteger(IL_IMAGE_WIDTH);

    // Get a pointer to the color component data. Although all rows of the image are accessible,
    // we'll need only the first row. The rest will be ignored
    ILubyte* imageData = ilGetData();

    if (imageFormat == IL_RGB) {
        // If there is no alpha channel in the file, we have to add an alpha value of 255 at every
        // fourth position. In order to do this, we'll go through both arrays at different paces
        // and add the required 255 in the 'newData' array
        ILubyte* newImageData = new ILubyte[dimensions_.x*4];
        for (int i = 0; i < dimensions_.x; ++i) {
            newImageData[i*4 + 0] = imageData[i*3 + 0];
            newImageData[i*4 + 1] = imageData[i*3 + 1];
            newImageData[i*4 + 2] = imageData[i*3 + 2];
            newImageData[i*4 + 3] = (char)(255);
        }
        imageData = newImageData;
    }
    LINFO("opened transfer function: " << filename);

    // Generate the keys from the given data
    generateKeys(reinterpret_cast<unsigned char*>(imageData));

    // Little cleanup
    ilDeleteImages(1, &imageID);

    // If we got up to here, everything was fine
    return true;
#else
bool TransFunc1DKeys::loadImage(const std::string&) {
    LERROR("Loading failed: No DevIL support.");
    return false;
#endif // VRN_MODULE_DEVIL
}

bool TransFunc1DKeys::loadTfi(const std::string& filename) {
    // open file stream
    std::ifstream stream(filename.c_str(), std::ios_base::in);
    if (stream.fail()) {
        LWARNING("Unable to open file " << filename << " for reading.");
        return false;
    }

    // deserialize from stream
    bool success = true;
    try {
        XmlDeserializer d(filename);
        d.read(stream);
        d.deserialize("TransFuncIntensity", *this);
        stream.close();
    }
    catch (SerializationException &e) {
        LWARNING("SerializationException: " << e.what());
        stream.close();
        success = false;
    }

    // log result
    if (success)
        LDEBUG("Loaded transfer function from file: " << filename);
    else
        LWARNING("Loading transfer function failed.");

    invalidateTexture();
    return success;
}

void TransFunc1DKeys::updateFrom(const TransFunc1DKeys& tf) {
    threshold_ = tf.threshold_;

    domain_ = tf.domain_;

    clearKeys();
    for (size_t i = 0; i < tf.keys_.size(); ++i) {
        TransFuncMappingKey* k = new TransFuncMappingKey(*(tf.keys_.at(i)));
        addKey(k);
    }
    gammaValue_ = tf.gammaValue_;
    alphaMode_ = tf.alphaMode_;
}

void TransFunc1DKeys::serialize(XmlSerializer& s) const {
    TransFunc::serialize(s);
    // serialize domain/threshold
    s.serialize("domain", domain_);
    s.serialize("threshold", threshold_);

    // serialize keys...
    s.serialize("Keys", keys_, "key");
}

void TransFunc1DKeys::deserialize(XmlDeserializer& s) {
    TransFunc::deserialize(s);

    // deserialize thresholds...
    s.optionalDeserialize("domain", domain_, tgt::vec2(0.0f, 1.0f));
    s.optionalDeserialize("threshold", threshold_,tgt::vec2(0.0f, 1.0f));

    // deserialize keys...
    s.deserialize("Keys", keys_, "key");
}

bool TransFunc1DKeys::loadTextTable(const std::string& filename) {
    if (filename == "")
        return false;

    FILE *file = fopen(filename.c_str(), "rt");

    if (file) {
        unsigned char* data = new unsigned char[256*4];
        for (int i = 0; i < 256; ++i) {
            int _red = 0, _green = 0, _blue = 0, _alpha = 0;
            int result = fscanf(file, "%d %d %d %d", &_red, &_green, &_blue, &_alpha);
            if (result != 4)
                LWARNING("loadTextTable(): invalid line");
            data[4*i + 0] = _red;
            data[4*i + 1] = _green;
            data[4*i + 2] = _blue;
            data[4*i + 3] = _alpha;
        }
        fclose(file);
        LINFO("opened transfer function: " << filename);
        dimensions_ = tgt::ivec3(256, 1, 1);
        generateKeys(data);
        delete[] data;

        return true;
    }
    else {
        LWARNING("unable to open transfer function: " << filename);
        return false;
    }
}

bool TransFunc1DKeys::loadOsirixCLUT(const std::string& filename) {
    LINFO("Opening Osirix CLUT: " << filename);

    TiXmlDocument doc(filename.c_str());

    //TODO: this loader is much to specific to a certain order of XML tags and will crash if
    // anything is missing (#276).
    if (doc.LoadFile()) {
        // read and check version of plist file
        TiXmlNode* currNode = doc.FirstChild("plist");
        TiXmlElement* currElement = currNode->ToElement();

        currNode = currNode->FirstChild("dict");
        currNode = currNode->FirstChild("key");
        currElement = currNode->ToElement();

        // get reference to red, green and blue channel
        TiXmlElement* blueElement = 0;
        TiXmlElement* greenElement = 0;
        TiXmlElement* redElement = 0;
        TiXmlNode* blueNode = currElement->NextSibling();
        TiXmlNode* greenNode = ((blueNode->NextSibling())->NextSibling());
        TiXmlNode* redNode = ((greenNode->NextSibling())->NextSibling());
        blueNode = blueNode->FirstChild("integer");
        greenNode = greenNode->FirstChild("integer");
        redNode = redNode->FirstChild("integer");

        unsigned char* data = new unsigned char[256*4];

        for (int i = 0; i < 256; ++i) {
            data[4*i + 0] = 0;
            data[4*i + 1] = 0;
            data[4*i + 2] = 0;
            data[4*i + 3] = 0;

            blueNode = blueNode->NextSibling("integer");
            greenNode = greenNode->NextSibling("integer");
            redNode = redNode->NextSibling("integer");

            if (blueNode == 0 || greenNode == 0 || redNode == 0)
                continue;

            blueElement = blueNode->ToElement();
            greenElement = greenNode->ToElement();
            redElement = redNode->ToElement();

            data[4*i + 0] = atoi(redElement->GetText());
            data[4*i + 1] = atoi(greenElement->GetText());
            data[4*i + 2] = atoi(blueElement->GetText());
            data[4*i + 3] = (char)(255);
        }

        dimensions_ = tgt::ivec3(256, 1, 1);
        generateKeys(data);
        delete[] data;

        return true;
    }
    else
        return false;
}

bool TransFunc1DKeys::loadImageJ(const std::string& filename) {
    std::ifstream fileStream(filename.c_str(), std::ios::in|std::ios::binary);

    if (fileStream.good()) {
        // Determine the length of file; necessary for typeswitch later
        size_t length;
        fileStream.seekg(0, std::ios::end);
        length = static_cast<size_t>(fileStream.tellg());
        fileStream.seekg(std::ios::beg);

        // If the length is that high, we can't use the data anyway
        if (length > 10000) {
            LWARNING("unable to open transfer function: " << filename);
            return false;
        }

        // The following procedures return 0 if the loading attempt was without success
        int size = 0;
        if (length > 768)
            // load NIH Image LUT
            size = openImageJBinary(fileStream, false);
        if ((size == 0) && ((length == 0) || (length == 768) || (length == 970)))
            //load raw LUT
            size = openImageJBinary(fileStream, true);
        if ((size == 0) && (length >= 768))
            size = openImageJText(fileStream);
        if (size == 0) {
            LWARNING("unable to open transfer function: " << filename);
        }

        // if size = 256, at least one of the procedures above was successful
        return (size == 256);

    }
    else {
        // The file couldn't be opened
        LWARNING("unable to open transfer function: " << filename);
        return false;
    }
}

int TransFunc1DKeys::openImageJBinary(std::ifstream& fileStream, bool raw) {
    // the default quantity of colors
    int numColors = 256;
    if (!raw) {
        // read the header information seperatly
        int id = readInt(fileStream);
        if (id != 1229147980) {
            // a leading 1229147980 (= 'ICOL') indicates an NIH Image LUT
            // if ICOL isn't the first entry in the table, let the other procedures handle it
            // but first, go back to the beginning
            fileStream.seekg(std::ios::beg);
            return 0;
        }
        readShort(fileStream); // Version
        numColors = readShort(fileStream); // Number of Colors
        readShort(fileStream);  // Start
        readShort(fileStream);  // End
        readDouble(fileStream); // Filler1
        readDouble(fileStream); // Filler2
        readInt(fileStream);    // Filler3
    }

    // The colors in a binary table are saved in succession so
    // first load the reds, then greens and at last blues
    char redColors[256];
    char greenColors[256];
    char blueColors[256];
    try {
        fileStream.read(&redColors[0], numColors);
        fileStream.read(&greenColors[0], numColors);
        fileStream.read(&blueColors[0], numColors);
    }
    catch (...) {
        throw;
    }

    unsigned char data[256*4];

    for (int i = 0; i < 256; ++i) {
        data[i*4 + 0] = redColors[i];
        data[i*4 + 1] = greenColors[i];
        data[i*4 + 2] = blueColors[i];
        data[i*4 + 3] = (char)(255);
    }

    dimensions_ = tgt::ivec3(256, 1, 1);
    generateKeys(&data[0]);

    return 256;
}

int TransFunc1DKeys::openImageJText(std::ifstream& fileStream) {
    // in this array, the converted values are stored
    // The maximum are 256*4 entries
    unsigned char* data = new unsigned char[256*4];
    int tmp;
    int numValues = 0;
    int strToIntConversionStatus = 0;

    while (strToIntConversionStatus != -1) {
        char* entry = new char[50];
        try {
            fileStream >> entry;

            // After this line StrToIntConversionStatus is:
            // 1 if 'line' was an integer ; 0 otherwise ; -1 if the end of the sfile is reached
            strToIntConversionStatus = sscanf(entry, "%i", &tmp);
        }
        catch (const std::exception& e) {
            LERROR(e.what());
            delete[] entry;
            return 0;
        }

        delete[] entry;

        // If an integer is encountered, add it to the array and increase the number of values
        if (strToIntConversionStatus == 1) {
            data[numValues] = tmp;
            ++numValues;
        }
    }

    /*  Now, two values of 'numValues' are possible. 256*3 = 768 or 256*4 = 1024.
        768, if there were 3 entries in a row
        1024, if there were 4 entires
    */
    if ((numValues != 768) && (numValues != 1024)) {
        // A wrong file was loaded
        delete[] data;
        return 0;
    }

    if (numValues == 1024) {
        // The first entry is an index, therefore it can be overwritten
        // Additional to that, an alpha value of 255 is inserted at the 4th position
        for (int i = 0; i < 256; ++i) {
            data[i*4 + 0] = data[i*4 + 1];
            data[i*4 + 1] = data[i*4 + 2];
            data[i*4 + 2] = data[i*4 + 3];
            data[i*4 + 3] = (char)(255);
        }

        dimensions_ = tgt::ivec3(256, 1, 1);
        generateKeys(data);
    }
    else {
        /* This one is a little bit tricky. We had 768 values in the file and added them
        sequentially to the data-array. Now we have to insert a '255' at every fourth position.
        In order to do that, we create a new array and move through both arrays at different
        paces.
        */
        unsigned char* newData = new unsigned char[1024];

        for (int i = 0; i < 256; ++i) {
            newData[i*4 + 0] = data[i*3 + 0];
            newData[i*4 + 1] = data[i*3 + 1];
            newData[i*4 + 2] = data[i*3 + 2];
            newData[i*4 + 3] = (char)(255);
        }

        dimensions_ = tgt::ivec3(256, 1, 1);
        generateKeys(newData);
        delete[] newData;
    }

    delete[] data;
    return 256;
}

void TransFunc1DKeys::generateKeys(unsigned char* data) {
    /* A short overview about the idea behind this method. For the sake of simplicity this is
    demonstrated with only one color channel, say 'Red'. It is generalized in the code below:

    We want to detect the peaks(=extrema) in the graph containing all the 'red'-values.
    In order to do this, we look at one 1/width-th at a time and compare the difference between
    the (i-2)th and (i-1)th point (= oldDelta_x) to the difference between the (i-1)th and ith point
    (stored in newDelta_x). Several possible things can happen:
    i) The difference doesn't change at all. This means that the difference between the points is
       linearly dependent and in this case, we don't have to add another key because we get linear
       interpolation from the the methods
    ii) The difference might be not-zero. This means the graph is discontinuous at this point and
        we have to insert a splitted point. On the left side we take the color from the
        (i-1)th point and on the right side from the ith point. The discontinuity is represented
        by the "jump" between the two parts of the splitted mapping key.
    iii) The difference could have changed by a factor of -1. In this case, we have to insert a
         mapping key right at the peak (i.e. at the ith point).

    If a key is placed, we don't want to place a key at the next location, because the difference
    will change not matter if there is a peak or not. There are a lot of cases, in this
    redundant keys will be generated and we don't want that


    In the code below, we generate a key whenever any of the colorchannels meets a criterion above.
    */

    // Storage for the old values
    int oldDeltaRed;
    int oldDeltaGreen;
    int oldDeltaBlue;
    int oldDeltaAlpha;

    // Storage for the new values
    int newDeltaRed;
    int newDeltaGreen;
    int newDeltaBlue;
    int newDeltaAlpha;

    // We want at least 2 values in the data array
    if (dimensions_.x < 2)
        return;

    clearKeys();

    addKey(new TransFuncMappingKey(0.f,
           col4(data[0], data[1], data[2], data[3])));
    addKey(new TransFuncMappingKey(1.f,
           col4(data[4*(dimensions_.x-1)+0], data[4*(dimensions_.x-1)+1], data[4*(dimensions_.x-1)+2], data[4*(dimensions_.x-1)+3])));

    // Calculate the starting point
    newDeltaRed   = data[4*1 + 0] - data[4*0 + 0];
    newDeltaGreen = data[4*1 + 1] - data[4*0 + 1];
    newDeltaBlue  = data[4*1 + 2] - data[4*0 + 2];
    newDeltaAlpha = data[4*1 + 3] - data[4*0 + 3];

    // The main loop. We start at 2 because the value for 1 already has been calculated.
    for (int iter = 2; iter < dimensions_.x; ++iter) {
        // Backup the old values and generate the new ones.
        oldDeltaRed = newDeltaRed;
        oldDeltaGreen = newDeltaGreen;
        oldDeltaBlue = newDeltaBlue;
        oldDeltaAlpha = newDeltaAlpha;

        newDeltaRed   = data[4*iter + 0] - data[4*(iter-1) + 0];
        newDeltaGreen = data[4*iter + 1] - data[4*(iter-1) + 1];
        newDeltaBlue  = data[4*iter + 2] - data[4*(iter-1) + 2];
        newDeltaAlpha = data[4*iter + 3] - data[4*(iter-1) + 3];

        // Has the difference quotient changed in any color channel?
        bool differenceQuotientChanged = (
            (oldDeltaRed   != newDeltaRed)   ||
            (oldDeltaGreen != newDeltaGreen) ||
            (oldDeltaBlue  != newDeltaBlue)  ||
            (oldDeltaAlpha != newDeltaAlpha));

        // Is the difference quotient different from zero in any channel?
        bool differenceQuotientNotZero = (
            (newDeltaRed   != 0) ||
            (newDeltaGreen != 0) ||
            (newDeltaBlue  != 0) ||
            (newDeltaAlpha != 0));

        // Has the difference quotient tilted in all channel's?
        // Mind the & instead of |
        bool differenceQuotientTilted = (
            (oldDeltaRed   == -newDeltaRed)   &&
            (oldDeltaGreen == -newDeltaGreen) &&
            (oldDeltaBlue  == -newDeltaBlue)  &&
            (oldDeltaAlpha == -newDeltaAlpha));

        if (differenceQuotientChanged) {
            if (differenceQuotientNotZero) {
                // We want to put a splitted key here (see ii above)
                TransFuncMappingKey* newkey = new TransFuncMappingKey(iter/static_cast<float>(dimensions_.x-1) ,
                    col4( data[4*(iter-1) + 0], data[4*(iter-1) + 1], data[4*(iter-1) + 2], data[4*(iter-1) + 3] )
                    );
                newkey->setSplit(true);
                newkey->setColorR(tgt::col4(data[4*iter + 0], data[4*iter + 1], data[4*iter + 2], data[4*iter + 3]));
                addKey(newkey);
            }
            else if (differenceQuotientTilted) {
                // We want a single key at i-1 here (see iii above)
                addKey(
                    new TransFuncMappingKey((iter - 1)/static_cast<float>(dimensions_.x-1),
                    col4(data[4*(iter-1) + 0], data[4*(iter-1) + 1], data[4*(iter-1) + 2], data[4*(iter-1) + 3])
                    ));
            }
            else {
                // Just add a key
                addKey(
                    new TransFuncMappingKey(iter/static_cast<float>(dimensions_.x-1),
                    col4(data[4*iter + 0], data[4*iter + 1], data[4*iter + 2], data[4*iter + 3])));
            }
        }
    }
}

tgt::vec2 TransFunc1DKeys::getDomain(int dimension) const {
    tgtAssert(dimension == 0, "Only one dimension!");

    if(dimension == 0)
        return domain_;
    else
        return tgt::vec2(0.0f);
}

void TransFunc1DKeys::setDomain(tgt::vec2 domain, int dimension) {
    tgtAssert(dimension == 0, "Only one dimension!");

    if(dimension == 0)
        domain_ = domain;
}

const PreIntegrationTable* TransFunc1DKeys::getPreIntegrationTable(float samplingStepSize, size_t dimension, bool useIntegral, bool computeOnGPU) {
    if (computeOnGPU && !preIntegrationProgram_)
        preIntegrationProgram_ = ShdrMgr.loadSeparate("passthrough.vert", "preintegration/preintegration.frag", generateStandardShaderHeader(0), false);

    return preIntegrationTableMap_.getPreIntegrationTable(samplingStepSize,dimension,useIntegral,computeOnGPU, preIntegrationProgram_);
}

void TransFunc1DKeys::deleteTexture() {
    TransFunc::deleteTexture();
    preIntegrationTableMap_.clear();
}

void TransFunc1DKeys::invertKeys() {
    std::vector<TransFuncMappingKey*> oldkeys(keys_);
    keys_.clear();

    for(std::vector<TransFuncMappingKey*>::reverse_iterator it = oldkeys.rbegin(); it != oldkeys.rend(); it++) {
        TransFuncMappingKey* key = *it;
        key->setIntensity(1.f - key->getIntensity());
        float tmp;
        tmp = key->getAlphaL();
        key->setAlphaL(key->getAlphaR());
        key->setAlphaR(tmp);
        tgt::col4 tmp2 = key->getColorL();
        key->setColorL(key->getColorR());
        key->setColorR(tmp2);
        keys_.push_back(key);
    }
    invalidateTexture();
}

// --------------------------------------------------------------------------------
// some stream helper functions

/**
 * Used not only to extract neccessary information but also to fast-forward through the header
 * Buffersize is given static, because the procedure shouldn't be edited if 64bit compiler is used
 * (and therefor the sizes i.e. int, double might change)
 * These functions are 'inline', because it might be faster if they are inlined by the compiler
 * (okay, not quite _that_ convincing i admit).
 *
 * @param stream the already opened file stream used to extract the data
 * @return the read int value
 */
inline int readInt(std::ifstream& stream) {
    char* buffer = new char[4];

    // The bytes are inserted backwards because of a discrepancy of
    // most significant bit <-> least significant bit between Java and C++
    for (int i = 3; i >= 0; --i)
        stream >> buffer[i];

    int ret = *(reinterpret_cast<int*>(buffer));
    delete[] buffer;
    return ret;
}

/**
 * \sa TransFuncIntensity::readInt(std::ifstream& stream)
 */
inline short readShort(std::ifstream& stream) {
    char* buffer = new char[2];

    for (int i = 1; i >= 0; --i)
        stream >> buffer[i];

    short ret = *(reinterpret_cast<short*>(buffer));
    delete[] buffer;
    return ret;
}

/**
 * \sa TransFuncIntensity::readInt(std::ifstream& stream)
 */
inline double readDouble(std::ifstream& stream) {
    char* buffer = new char[8];

    for (int i = 7; i >= 0; --i)
        stream >> buffer[i];

    double ret = *(reinterpret_cast<double*>(buffer));
    delete[] buffer;
    return ret;
}

} // namespace voreen
