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

#include "octreecreator.h"

#include "voreen/core/datastructures/octree/volumeoctree.h"
#include "voreen/core/datastructures/octree/octreebrickpoolmanager.h"
#include "voreen/core/datastructures/octree/octreebrickpoolmanagerdisk.h"
#include "voreen/core/datastructures/octree/octreeutils.h"

#include "voreen/core/voreenapplication.h"
#include "voreen/core/datastructures/volume/volumeram.h"
#include "voreen/core/datastructures/volume/volumedisk.h"
#include "voreen/core/datastructures/volume/volume.h"
#include "voreen/core/datastructures/volume/operators/volumeoperatorresizepoweroftwo.h"
#include "voreen/core/datastructures/volume/operators/volumeoperatorconvert.h"
#include "voreen/core/datastructures/volume/volumedecorator.h"
#include "voreen/core/datastructures/meta/realworldmappingmetadata.h"
#include "voreen/core/datastructures/volume/volumeminmax.h"
#include "voreen/core/datastructures/volume/histogram.h"
#include "voreen/core/datastructures/volume/volumehash.h"
#include "voreen/core/datastructures/volume/volumepreview.h"
#include "voreen/core/datastructures/octree/octreebrickpoolmanagerdisk.h"
#include "voreen/core/datastructures/octree/octreeutils.h"
#include "voreen/core/io/serialization/serialization.h"
#include "voreen/core/utils/hashing.h"
#include "voreen/core/utils/memoryinfo.h"

#include "tgt/vector.h"
#include "tgt/tgt_math.h"
#include "tgt/stopwatch.h"
#include "tgt/filesystem.h"

#include <queue>
#include <stack>

using tgt::ivec2;
using tgt::ivec3;
using tgt::ivec4;
using tgt::vec2;
using tgt::vec3;
using tgt::vec4;
using tgt::svec3;

namespace {

const std::string CACHE_SUBDIR =             "OctreeCreator";
const std::string OCTREE_FILENAME =          "octree.xml";
const std::string BRICK_BUFFER_SUBDIR =      "brickBuffer";
const std::string BRICK_BUFFER_FILE_PREFIX = "buffer_";

/// Retrieves the octree cache limit from the VoreenApplication.
uint64_t getCacheLimit() {
    const std::string propertyID = "octreeCacheLimit";
    tgtAssert(voreen::VoreenApplication::app(), "VoreenApplication not instantiated");
    voreen::IntProperty* cacheSizeProp = dynamic_cast<voreen::IntProperty*>(voreen::VoreenApplication::app()->getProperty(propertyID));
    if (cacheSizeProp) {
        return (uint64_t)cacheSizeProp->get() * (1<<30);
    }
    else {
        LWARNINGC("voreen.OctreeCreator", "IntProperty '" << propertyID << "' not found");
        return 0;
    }
}

/// Returns the brick dimension for a given volume dimension and tree depth.
size_t computeBrickDim(const tgt::svec3& volumeDim, size_t numLevels) {
    tgtAssert(numLevels > 0, "numLevels is 0");
    size_t octreeDim = tgt::nextLargerPowerOfTwo(tgt::max(volumeDim));
    size_t brickDim = octreeDim / (size_t)(1 << (numLevels-1));
    return brickDim;
}

/// Returns the octree depth for the passed brick and volume dimensions
size_t computeTreeDepth(const tgt::svec3& brickDim, const tgt::svec3& volumeDim) {
    tgtAssert(isCubicAndPot(brickDim), "brick dim not cubic-pot");
    size_t octreeDim = tgt::nextLargerPowerOfTwo(tgt::max(volumeDim));
    size_t treeDepth = tgt::ilog2((int)octreeDim / (int)brickDim.x) + 1;
    return treeDepth;
}

size_t computeMaxTreeDepth(const tgt::svec3& volumeDim) {
    size_t octreeDim = tgt::nextLargerPowerOfTwo(tgt::max(volumeDim));
    return (size_t)tgt::ilog2((int)octreeDim);
}

} // namespace anonymous

namespace voreen {

const std::string OctreeCreator::loggerCat_("voreen.OctreeCreator");

OctreeCreator::OctreeCreator()
    : VolumeProcessor()
    , volumeInport_(Port::INPORT, "volumeInput", "Volume Input")
    , volumeInport2_(Port::INPORT, "volumeInput2", "Volume Input 2")
    , volumeInport3_(Port::INPORT, "volumeInput3", "Volume Input 3")
    , volumeInport4_(Port::INPORT, "volumeInput4", "Volume Input 4")
    , volumeOutport_(Port::OUTPORT, "volumeOutport", "Volume Output")
    , generateOctreeButton_("generateOctree", "Generate/Load Octree")
    , statusProperty_("statusProperty", "", "", VALID)
    , progressProperty_("progressProperty", "")
    , autogenerateOctree_("autogenerateOctree", "Auto-generate", false)
    , brickDimensions_("brickDimensions", "Brick Dimensions")
    , treeDepth_("treeDepth", "Tree Depth", 6, 1, 8)
    , homogeneityThreshold_("homogeneityThreshold", "Homogeneity Threshold", 0.05f, 0.f, 0.1f)
    , useRelativeThreshold_("useRelativeThreshold", "Relative Threshold", true)
    , brickPoolManager_("brickPoolManager", "Brick Pool Manager")
    , singleBufferMemorySize_("singleBufferMemorySize", "Page File Size (MB)", 32, 1, 256)
    , numThreads_("numThreads", "Num Threads", 8, 1, 16, VALID)
    , clearOctree_("clearOctree", "Clear Octree")
    , forceGenerate_(false)
    , currentConfigurationHash_("")
{
    addPort(volumeInport_);
    addPort(volumeInport2_);
    addPort(volumeInport3_);
    addPort(volumeInport4_);
    addPort(volumeOutport_);

    statusProperty_.setReadOnly(true);
    addProperty(generateOctreeButton_);
    generateOctreeButton_.onChange(CallMemberAction<OctreeCreator>(this, &OctreeCreator::forceRegenerate));
    addProperty(statusProperty_);
    addProperty(progressProperty_);
    addProgressBar(&progressProperty_);
    addProperty(autogenerateOctree_);
    autogenerateOctree_.onChange(CallMemberAction<OctreeCreator>(this, &OctreeCreator::updatePropertyConfiguration));

    brickDimensions_.addOption("treeDepth", "Derive from Tree Depth",   0);
    //brickDimensions_.addOption("2",       "[2 2 2]",                  2);
    brickDimensions_.addOption("4",         "[4 4 4]",                  4);
    brickDimensions_.addOption("8",         "[8 8 8]",                  8);
    brickDimensions_.addOption("16",        "[16 16 16]",               16);
    brickDimensions_.addOption("32",        "[32 32 32]",               32);
    brickDimensions_.addOption("64",        "[64 64 64]",               64);
    brickDimensions_.addOption("128",       "[128 128 128]",            128);
    brickDimensions_.selectByKey("treeDepth");
    brickDimensions_.onChange(CallMemberAction<OctreeCreator>(this, &OctreeCreator::updatePropertyConfiguration));
    addProperty(brickDimensions_);
    treeDepth_.setWidgetsEnabled(false);
    addProperty(treeDepth_);
    //homogeneityThreshold_.setNumDecimals(2);
    //homogeneityThreshold_.setStepping(0.01f);
    homogeneityThreshold_.setTracking(false);
    addProperty(homogeneityThreshold_);
    addProperty(useRelativeThreshold_);
    treeDepth_.setGroupID("configuration");
    brickDimensions_.setGroupID("configuration");
    homogeneityThreshold_.setGroupID("configuration");
    useRelativeThreshold_.setGroupID("configuration");

    brickPoolManager_.addOption("brickPoolManagerRAM",  "RAM (non-persistent)");
    brickPoolManager_.addOption("brickPoolManagerDisk", "Disk");
    brickPoolManager_.select("brickPoolManagerDisk");
    addProperty(brickPoolManager_);
    brickPoolManager_.onChange(CallMemberAction<OctreeCreator>(this, &OctreeCreator::updatePropertyConfiguration));
    addProperty(singleBufferMemorySize_);
    brickPoolManager_.setGroupID("configuration");
    singleBufferMemorySize_.setGroupID("configuration");
    numThreads_.setGroupID("configuration");
    addProperty(numThreads_);
    setPropertyGroupGuiName("configuration", "Octree Configuration");

    addProperty(clearOctree_);
    clearOctree_.onChange(CallMemberAction<OctreeCreator>(this, &OctreeCreator::clearOctree));
}

OctreeCreator::~OctreeCreator() {
}

Processor* OctreeCreator::create() const {
    return new OctreeCreator();
}

void OctreeCreator::initialize() throw (VoreenException) {
    VolumeProcessor::initialize();

    updatePropertyConfiguration();

    limitCacheSize(VoreenApplication::app()->getCachePath(), getCacheLimit(), false);
}

void OctreeCreator::deinitialize() throw (VoreenException) {
    volumeOutport_.clear();
    if (!VoreenApplication::app()->useCaching() && currentConfigurationHash_ != "")
        clearOctree();

    limitCacheSize(VoreenApplication::app()->getCachePath(), getCacheLimit(), false);

    VolumeProcessor::deinitialize();
}

void OctreeCreator::adjustPropertiesToInput() {
    updatePropertyConfiguration();
}

bool OctreeCreator::isReady() const {
    return volumeInport_.isReady() && volumeOutport_.isReady();
}

void OctreeCreator::process() {
    if (!volumeOutport_.hasData())
        currentConfigurationHash_ = "";

    std::string configHash = getConfigurationHash();

    if (configHash == currentConfigurationHash_ && !forceGenerate_)
        updatePropertyConfiguration();

    volumeOutport_.clear();
    if (!VoreenApplication::app()->useCaching() && currentConfigurationHash_ != "")
        clearOctree();
    currentConfigurationHash_ = "";

    VolumeBase* inputVolume = const_cast<VolumeBase*>(volumeInport_.getData());
    if (!inputVolume || (!autogenerateOctree_.get() && !forceGenerate_)) {
        updatePropertyConfiguration();
        return;
    }

    // remove existing octree representation
    inputVolume->removeRepresentation<VolumeOctreeBase>();

    // if caching enabled, try to restore octree from cache
    VolumeOctreeBase* octree = 0;
    if (VoreenApplication::app()->useCaching()) {
        octree = restoreOctreeFromCache();
        if (octree)
            statusProperty_.set("Octree loaded from cache.");
    }

    // if octree has not been restored from cache, generate it and store it to cache
    if (!octree && (autogenerateOctree_.get() || forceGenerate_)) {
        statusProperty_.set("Generating octree...");
        tgt::Stopwatch watch;
        watch.start();
        try {
            octree = generateOctree();
        }
        catch (VoreenException& e) {
            std::string errorMsg = "Failed to generate octree: " + std::string(e.what());
            LERROR(errorMsg);
            VoreenApplication::app()->showMessageBox(getGuiName(), errorMsg, true);
            statusProperty_.set("Failed to generate octree!");
            forceGenerate_ = false;
            return;
        }
        tgtAssert(octree, "no octree"); //< otherwise exception expected

        statusProperty_.set("Octree generated.");
        LINFO("Octree construction time: " << formatTime(watch.getRuntime()));
        LDEBUG("- " << MemoryInfo::getProcessMemoryUsageAsString());
        LDEBUG("- " << MemoryInfo::getAvailableMemoryAsString());

        // store octree to cache, if enabled
        if (octree && VoreenApplication::app()->useCaching() && !brickPoolManager_.isSelected("brickPoolManagerRAM")) {
            storeOctreeToCache(octree);
        }
    }
    forceGenerate_ = false;

    if (octree) {
        currentConfigurationHash_ = configHash;

        // assign RAM limit
        size_t ramLimit = VoreenApplication::app()->getCpuRamLimit();
        if (const OctreeBrickPoolManagerDisk* brickPoolManager =
            dynamic_cast<const OctreeBrickPoolManagerDisk*>(static_cast<VolumeOctree*>(octree)->getBrickPoolManager())) {
                const_cast<OctreeBrickPoolManagerDisk*>(brickPoolManager)->setRAMLimit(ramLimit);
        }

        octree->logDescription();

        // min/max values
        std::vector<float> minValues, maxValues, minNormValues, maxNormValues;
        for (size_t i=0; i<octree->getNumChannels(); i++) {
            float minNorm = octree->getRootNode()->getMinValue(i) / 65535.f;
            float maxNorm = octree->getRootNode()->getMaxValue(i) / 65535.f;
            tgtAssert(minNorm <= maxNorm, "invalid min/max values");
            float min = inputVolume->getRealWorldMapping().normalizedToRealWorld(minNorm);
            float max = inputVolume->getRealWorldMapping().normalizedToRealWorld(maxNorm);
            minValues.push_back(min);
            maxValues.push_back(max);
            minNormValues.push_back(minNorm);
            maxNormValues.push_back(maxNorm);
        }
        VolumeMinMax* volumeMinMax = new VolumeMinMax(minValues, maxValues, minNormValues, maxNormValues);

        // histograms
        std::vector<Histogram1D> histograms;
        for (size_t i=0; i<octree->getNumChannels(); i++) {
            histograms.push_back(Histogram1D(*(octree->getHistogram(i))));
        }
        VolumeHistogramIntensity* histogramData = new VolumeHistogramIntensity(histograms);

        // single-channel uint16_t octree => add octree to input volume handle
        if (octree->getNumChannels() == 1 && inputVolume->getFormat() == "uint16") {
            // assign minmax derived data to inputVolume volume, if not present
            if (!inputVolume->hasDerivedData<VolumeMinMax>())
                const_cast<VolumeBase*>(inputVolume)->addDerivedData(volumeMinMax);
            else
                delete volumeMinMax;

            // assign histograms, if not present
            if (!inputVolume->hasDerivedData<VolumeHistogramIntensity>())
                const_cast<VolumeBase*>(inputVolume)->addDerivedData(histogramData);
            else
                delete histogramData;

            inputVolume->addRepresentation(octree);

            volumeOutport_.setData(inputVolume, false);
        }
        else { // multi-channel volume or differing data type => create new volume handle
            VolumeBase* outputVolume = new Volume(octree, inputVolume);

            outputVolume->addDerivedData(volumeMinMax);
            outputVolume->addDerivedData(histogramData);
            outputVolume->getDerivedData<VolumePreview>(); //< prevent creation in background thread (brick pool access)

            // TODO: add hash computation to VolumeOctreeBase interface
            VolumeHash* volumeHash = new VolumeHash();
            volumeHash->setHash(VoreenHash::getHash(getConfigurationHash()));
            outputVolume->addDerivedData(volumeHash);

            volumeOutport_.setData(outputVolume, true);
        }
        tgtAssert(volumeOutport_.hasData(), "outport has no data");

        setProgress(1.f);
    }
    else {
        setProgress(0.f);
    }

    updatePropertyConfiguration();
}

VolumeOctreeBase* OctreeCreator::generateOctree() throw (VoreenException) {
    // retrieve input RAM volumes
    std::vector<const VolumeBase*> inputVolumes;
    if (volumeInport_.hasData())
        inputVolumes.push_back(volumeInport_.getData());
    if (volumeInport2_.hasData())
        inputVolumes.push_back(volumeInport2_.getData());
    if (volumeInport3_.hasData())
        inputVolumes.push_back(volumeInport3_.getData());
    if (volumeInport4_.hasData())
        inputVolumes.push_back(volumeInport4_.getData());

    if (inputVolumes.empty())
        throw VoreenException("No input volumes");

    const tgt::svec3 volumeDim = inputVolumes.front()->getDimensions();

    // select brick pool manager
    OctreeBrickPoolManagerBase* brickPoolManager = 0;
    if (brickPoolManager_.isSelected("brickPoolManagerRAM"))
        brickPoolManager = new OctreeBrickPoolManagerRAM(singleBufferMemorySize_.get() << 20);
    else if (brickPoolManager_.isSelected("brickPoolManagerDisk")) {
        std::string brickPoolPath = tgt::FileSystem::cleanupPath(getOctreeStoragePath() + "/" + BRICK_BUFFER_SUBDIR);
        if (!tgt::FileSystem::dirExists(brickPoolPath))
            tgt::FileSystem::createDirectoryRecursive(brickPoolPath);
        brickPoolManager = new OctreeBrickPoolManagerDisk(static_cast<size_t>(singleBufferMemorySize_.get()) << 20,
            VoreenApplication::app()->getCpuRamLimit(), brickPoolPath, BRICK_BUFFER_FILE_PREFIX);
    }
    else {
        throw VoreenException("Unknown brick pool manager: " + brickPoolManager_.get());
    }

    // determine homogeneity threshold
    float homogeneityThreshold = homogeneityThreshold_.get();
    if (useRelativeThreshold_.get()) {
        // use middle slice of first volume to calibrate threshold
        size_t sliceID = volumeDim.z / 2;
        try {
            VolumeRAM* slice = 0;
            if (inputVolumes.front()->hasRepresentation<VolumeRAM>()) {
                slice = inputVolumes.front()->getRepresentation<VolumeRAM>()->getSubVolume(tgt::svec3(volumeDim.x, volumeDim.y, 1), tgt::svec3(0, 0, sliceID));
            }
            else if (inputVolumes.front()->hasRepresentation<VolumeDisk>()) {
                slice = inputVolumes.front()->getRepresentation<VolumeDisk>()->loadSlices(sliceID, sliceID);
            }
            tgtAssert(slice, "no slice");

            float min = slice->minNormalizedValue();
            float max = slice->maxNormalizedValue();
            tgtAssert(min <= max, "invalid min/max values");
            homogeneityThreshold *= (max-min);

            delete slice;

            LDEBUG("Middle slice min/max: " << min << "/" << max);
            LDEBUG("Homogeneity threshold: " << homogeneityThreshold);
        }
        catch (std::exception& e) {
            LWARNING("Failed to extract slice " << sliceID << " from input volume for min/max range computation: " << e.what());
        }
    }

    // create octree
    VolumeOctreeBase* octree = 0;
    try {
        size_t brickDim;
        if (brickDimensions_.isSelected("treeDepth")) {
            brickDim = computeBrickDim(inputVolumes.front()->getDimensions(), treeDepth_.get());
            brickDim = std::max<size_t>(brickDim, 4);
        }
        else {
            brickDim = brickDimensions_.getValue();
        }
        tgtAssert(brickDim >= 4, "brick dim too small");
        tgtAssert(tgt::isPowerOfTwo((int)brickDim), "brick dim not power-of-two");

        octree = new VolumeOctree(inputVolumes, (int)brickDim, homogeneityThreshold, brickPoolManager, numThreads_.get(), this);
    }
    catch (std::exception& e) {
        LINFO(MemoryInfo::getProcessMemoryUsageAsString());
        LINFO(MemoryInfo::getAvailableMemoryAsString());
        throw VoreenException(e.what());
    }

    return octree;
}

void OctreeCreator::storeOctreeToCache(const VolumeOctreeBase* octree) const {
    tgtAssert(octree, "null pointer passed");

    // detect/create base cache path
    const std::string cachePath = getOctreeStoragePath();
    if (!tgt::FileSystem::dirExists(cachePath) && !tgt::FileSystem::createDirectoryRecursive(cachePath)) {
        LWARNING("Failed to create cache directory: " << cachePath);
        return;
    }
    else {
        //tgt::FileSystem::clearDirectory(cachePath);
    }

    // serialize octree
    std::string octreePath = tgt::FileSystem::cleanupPath(cachePath + "/" + OCTREE_FILENAME);
    XmlSerializer serializer(octreePath);
    try {
        serializer.serialize("Octree", octree);
    }
    catch (SerializationException& e) {
        LWARNING("Failed to serialize octree: " << e.what());
        return;
    }

    // write serialization to stream
    std::ostringstream textStream;
    try {
        serializer.write(textStream);
        if (textStream.fail()) {
            LWARNING("Failed to write octree serialization to string stream");
            return;
        }
    }
    catch (std::exception& e) {
        LWARNING("Failed to write octree serialization to string stream: " << e.what());
        return;
    }

    // now we have a valid string stream containing the serialized octree
    // => open output file and write it to the file
    std::fstream fileStream(octreePath.c_str(), std::ios_base::out);
    if (fileStream.fail()) {
        LWARNING("Failed to open file '" << octreePath << "' for writing.");
        return;
    }

    try {
        fileStream << textStream.str();
    }
    catch (std::exception& e) {
        LWARNING("Failed to write serialization data stream to file '" << octreePath << "': " << std::string(e.what()));
    }
    fileStream.close();

    LINFO("Serialized octree to file: " << octreePath);

    limitCacheSize(VoreenApplication::app()->getCachePath(), getCacheLimit(), true);
}

VolumeOctreeBase* OctreeCreator::restoreOctreeFromCache() const {
    const std::string octreeFile = tgt::FileSystem::cleanupPath(getOctreeStoragePath() + "/" + OCTREE_FILENAME);
    if (!tgt::FileSystem::fileExists(octreeFile)) {
        LINFO("No matching octree found in cache.");
        return 0;
    }

    tgt::Stopwatch stopWatch;
    stopWatch.start();

    // open file for reading
    std::fstream fileStream(octreeFile.c_str(), std::ios_base::in);
    if (fileStream.fail()) {
        LWARNING("Failed to open cached octree file '" << octreeFile << "' for reading.");
        return 0;
    }

    // read data stream into deserializer
    XmlDeserializer d(octreeFile);
    try {
        d.read(fileStream);
    }
    catch (std::exception& e) {
        LWARNING("Failed to read serialization data stream from cached octree file '" + octreeFile + "': " + e.what());
        return 0;
    }

    // deserialize octree from data stream
    VolumeOctreeBase* octree = 0;
    try {
        LDEBUG("Restoring cached octree from file: " << octreeFile );
        LDEBUG("- Before: " << MemoryInfo::getProcessMemoryUsageAsString());
        LDEBUG("- Before: " << MemoryInfo::getAvailableMemoryAsString());

        d.deserialize("Octree", octree);
        tgtAssert(octree, "null pointer after deserialization");
        tgtAssert(octree->getRootNode(), "deserialized octree has no root node");

        LINFO("Restored cached octree from file: " << octreeFile << " (" << stopWatch.getRuntime() << " msec)");
        LDEBUG("- After: " << MemoryInfo::getProcessMemoryUsageAsString());
        LDEBUG("- After: " << MemoryInfo::getAvailableMemoryAsString());

        // update file access time
        tgt::FileSystem::updateFileTime(octreeFile);
    }
    catch (std::exception& e) {
        LWARNING("Failed to restore cached octree from file '" + octreeFile + "': " + e.what());
        delete octree;
        octree = 0;
    }

    return octree;
}

void OctreeCreator::clearOctree() {
    // clear octree
    volumeOutport_.clear();
    currentConfigurationHash_ = "";

    // clear disk cache
    std::string octreePath = getOctreeStoragePath();
    if (tgt::FileSystem::dirExists(octreePath)) {
        if (!tgt::FileSystem::deleteDirectoryRecursive(octreePath))
            LWARNING("Failed to delete octree on disk: " + octreePath);
    }

    //setProgress(0.f);
    updatePropertyConfiguration();
}

std::string OctreeCreator::getOctreeStoragePath() const {
    return VoreenApplication::app()->getCachePath(CACHE_SUBDIR + "/" + getConfigurationHash());
}

std::string OctreeCreator::getConfigurationHash() const {
    // get input volume hash
    std::string volumeHash;
    if (volumeInport_.hasData())
        volumeHash += volumeInport_.getData()->getHash();
    if (volumeInport2_.hasData())
        volumeHash += volumeInport2_.getData()->getHash();
    if (volumeInport3_.hasData())
        volumeHash += volumeInport3_.getData()->getHash();
    if (volumeInport4_.hasData())
        volumeHash += volumeInport4_.getData()->getHash();
    volumeHash = VoreenHash::getHash(volumeHash);

    // compute property string
    std::map<std::string, const Property*> propertyMap;
    propertyMap[brickDimensions_.getID()] = &brickDimensions_;
    propertyMap[treeDepth_.getID()] = &treeDepth_;
    propertyMap[homogeneityThreshold_.getID()] = &homogeneityThreshold_;
    propertyMap[useRelativeThreshold_.getID()] = &useRelativeThreshold_;
    propertyMap[brickPoolManager_.getID()] = &brickPoolManager_;
    propertyMap[singleBufferMemorySize_.getID()] = &singleBufferMemorySize_;

    XmlSerializer s;
    const bool usePointerContentSerialization = s.getUsePointerContentSerialization();
    s.setUsePointerContentSerialization(true);
    try {
        s.serialize("Properties", propertyMap, "Property", "name");
    }
    catch (SerializationException& e) {
        LWARNING(e.what());
    }
    s.setUsePointerContentSerialization(usePointerContentSerialization);
    std::stringstream propertyConfig;
    s.write(propertyConfig);
    std::string propertyHash = VoreenHash::getHash(propertyConfig.str());

    // concatenate property and volume hash
    return volumeHash + "-" + propertyHash;
}

void OctreeCreator::forceRegenerate() {
    forceGenerate_ = true;
    invalidate();
}

void OctreeCreator::updatePropertyConfiguration() {
    //generateOctreeButton_.setWidgetsEnabled(!autogenerateOctree_.get());

    treeDepth_.setWidgetsEnabled(brickDimensions_.isSelected("treeDepth"));

    if (volumeInport_.hasData()) {
        tgt::svec3 volumeDim = volumeInport_.getData()->getDimensions();

        /*size_t maxDepth = tgt::clamp(computeMaxTreeDepth(volumeDim), (size_t)1, (size_t)8);
        treeDepth_.setMaxValue((int)maxDepth);

        size_t brickDim = computeBrickDim(volumeDim, treeDepth_.get());
        tgtAssert(tgt::isPowerOfTwo((int)brickDim), "brick dim not power-of-two");
        brickDimensions_.set((int)brickDim); */

        if (!brickDimensions_.isSelected("treeDepth")) {
            tgt::svec3 brickDim(brickDimensions_.getValue());
            size_t treeDepth = computeTreeDepth(brickDim, volumeDim);
            treeDepth_.set((int)treeDepth);
        }

        // is octree present in cache?
        bool inCache = false;
        if (VoreenApplication::app()->useCaching()) {
            std::string octreeFile = tgt::FileSystem::cleanupPath(getOctreeStoragePath() + "/octree.xml");
            inCache = tgt::FileSystem::fileExists(octreeFile);
        }

        // generate button
        generateOctreeButton_.setWidgetsEnabled(!volumeOutport_.hasData() || currentConfigurationHash_ != getConfigurationHash());
        generateOctreeButton_.setGuiName(inCache ? "Load Octree" : "Generate Octree");
        generateOctreeButton_.updateWidgets();

        // status property
        if (!volumeOutport_.hasData()) {
            progressProperty_.setProgress(0.f);
            statusProperty_.set(inCache ? "Octree found in cache." : "Octree not found in cache!");
        }

    }
    else {
        generateOctreeButton_.setWidgetsEnabled(false);
        progressProperty_.setProgress(0.f);
        statusProperty_.set("No input data.");
    }

}

// statics

void OctreeCreator::limitCacheSize(const std::string& cacheDirectory, const uint64_t maxCacheSize, bool keepLatest) {
    std::string octreeCachePath = tgt::FileSystem::cleanupPath(cacheDirectory + "/" + CACHE_SUBDIR);
    if (!tgt::FileSystem::dirExists(octreeCachePath))
        return;

    if (maxCacheSize == 0)
        return;

    uint64_t cacheSize = tgt::FileSystem::dirSize(octreeCachePath);
    if (cacheSize <= maxCacheSize)
        return;

    LINFO("Octree cache size (" << formatMemorySize(cacheSize) << ") exceeds limit (" << formatMemorySize(maxCacheSize) << "). " << "Cleaning cache...");

    // collect octree directories and their access time in chronologically ordered list
    std::list< std::pair<std::string, time_t> > octreeDirs;
    std::vector<std::string> subDirs = tgt::FileSystem::listSubDirectories(octreeCachePath);
    for (size_t i=0; i<subDirs.size(); i++) {
        std::string octreeDir = tgt::FileSystem::cleanupPath(octreeCachePath + "/" + subDirs.at(i));
        std::string octreeFile = octreeDir + "/" + OCTREE_FILENAME;
        if (tgt::FileSystem::fileExists(octreeFile)) {
            time_t octreeTime = tgt::FileSystem::fileAccessTime(octreeFile);
            for (std::list< std::pair<std::string, time_t> >::iterator it = octreeDirs.begin(); ;++it) {
                if (it == octreeDirs.end() || it->second > octreeTime) {
                    octreeDirs.insert(it, std::make_pair(octreeDir, octreeTime));
                    break;
                }
            }
        }
        else { // invalid directory (no octree.xml) => delete first
            octreeDirs.insert(octreeDirs.begin(), std::make_pair(octreeDir, 0));
        }
    }

    // delete collected octree dir in chronological order until cache does not exceed the limit anymore
    int minCacheElems = (keepLatest ? 1 : 0);
    while (cacheSize > maxCacheSize && (octreeDirs.size() > minCacheElems) ) {
        std::string octreeDir = octreeDirs.front().first;
        uint64_t dirSize = tgt::FileSystem::dirSize(octreeDir);
        LINFO("Deleting cached octree: " << octreeDir << " (size: " << formatMemorySize(dirSize) << ")");
        if (!tgt::FileSystem::deleteDirectoryRecursive(octreeDir))
            LWARNING("Failed to delete octree directory: " << octreeDir);
        octreeDirs.pop_front();
        cacheSize -= std::min(dirSize, cacheSize);
    }
}

void OctreeCreator::deleteCache(const std::string& cacheDirectory) {
    std::string octreeCachePath = tgt::FileSystem::cleanupPath(cacheDirectory + "/" + CACHE_SUBDIR);
    if (tgt::FileSystem::dirExists(octreeCachePath)) {
        LINFOC("voreen.OctreeCreator", "Clearing octree cache directory: " << octreeCachePath);
        if (!tgt::FileSystem::deleteDirectoryRecursive(octreeCachePath))
            LWARNINGC("voreen.OctreeCreator", "Failed to delete octree cache directory: " << octreeCachePath);
    }

}

}   // namespace
