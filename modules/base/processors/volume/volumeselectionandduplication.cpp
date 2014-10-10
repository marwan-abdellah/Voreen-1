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

#include "volumeselectionandduplication.h"

#include "voreen/core/datastructures/volume/volumeatomic.h"
#include "voreen/core/datastructures/volume/operators/volumeoperatorresize.h"
#include "voreen/core/datastructures/geometry/meshlistgeometry.h"

using tgt::ivec3;
using tgt::svec3;
using tgt::vec3;

namespace voreen {

const std::string VolumeSelectionAndDuplication::loggerCat_("voreen.VolumeSelectionAndDuplication");

VolumeSelectionAndDuplication::VolumeSelectionAndDuplication()
    : CachingVolumeProcessor()
    , inportFirst_(Port::INPORT, "volume.first", "Volume1 Input")
    , inportSecond_(Port::INPORT, "volume.second", "Volume2 Input")
    , outportFirst_(Port::OUTPORT, "outport.first", "Volume Output", true)
    , outportSecond_(Port::OUTPORT, "outport.second", "Volume Output", true)
    , enableProcessing_("enabled", "Enable", true)
    , selectedVolume_("selectedVolume", "Selected Volume")
    , filteringMode_("filteringMode", "Filtering")
{
    addPort(inportFirst_);
    addPort(inportSecond_);
    addPort(outportFirst_);
    addPort(outportSecond_);

    selectedVolume_.addOption("first",  "First",  FIRST_VOLUME);
    selectedVolume_.addOption("second", "Second", SECOND_VOLUME);
    selectedVolume_.select("first");

    filteringMode_.addOption("nearest", "Nearest");
    filteringMode_.addOption("linear",  "Linear ");
    filteringMode_.addOption("cubic",   "Cubic");
    filteringMode_.set("linear");

    addProperty(enableProcessing_);
    addProperty(filteringMode_);
    addProperty(selectedVolume_);
}

VolumeSelectionAndDuplication::~VolumeSelectionAndDuplication() {}

Processor* VolumeSelectionAndDuplication::create() const {
    return new VolumeSelectionAndDuplication();
}

void VolumeSelectionAndDuplication::process() {
    // Getting the input volumes.
    tgtAssert(inportFirst_.getData() && inportFirst_.getData()->getRepresentation<VolumeRAM>(), "No input volume");
    tgtAssert(inportSecond_.getData() && inportSecond_.getData()->getRepresentation<VolumeRAM>(), "No input volume");

    const VolumeBase* firstVolume = inportFirst_.getData();
    const VolumeBase* secondVolume = inportSecond_.getData();

    // The selected volume will be allocated and then the data of the selected one
    // will be added to it.
    Volume* selectedVolume = 0;

    // If not processingm just by-pass the first volume.
    if (!enableProcessing_.get()) {
        outportFirst_.setData(const_cast<VolumeBase*>(inportFirst_.getData()), false);
        outportSecond_.setData(const_cast<VolumeBase*>(inportFirst_.getData()), false);
        return;
    }

    // This is to accelerate the processing.
    if (firstVolume->getNumChannels() == secondVolume->getNumChannels()) {
        LINFO("Performing channel-wise combination.");
    }
    else {
        LINFO("Number of channels of input volumes do not match."
              "Combining each channel of second volume with channel 0 of first volume");
    }


    // optimized combination for volumes that share a common grid in world-space
    if (firstVolume->getDimensions() == secondVolume->getDimensions() &&
        firstVolume->getSpacing() == secondVolume->getSpacing()     &&
        firstVolume->getOffset() == secondVolume->getOffset()     &&
        firstVolume->getPhysicalToWorldMatrix() == secondVolume->getPhysicalToWorldMatrix()) {

        try {
            if (selectedVolume_.isSelected("first"))
                selectedVolume = firstVolume->clone();
            else
                selectedVolume = secondVolume->clone();
        }
        catch (const std::bad_alloc&) {
            LERROR("Failed to create combined volume with dimensions " << firstVolume->getDimensions()
                << " : bad allocation");
        }

        if (selectedVolume) {
            LINFO("Performing optimized combination on common grid with dimensions "
                << firstVolume->getDimensions() << "...");
            copyInputToSelectedVolumeOnCommonGrid(selectedVolume, firstVolume, secondVolume, selectedVolume_.getValue());
        }
    }
    // standard combination with resampling
    else {
        if (selectedVolume_.isSelected("first"))
            selectedVolume = createSelectedVolume(firstVolume, secondVolume);
        else
            selectedVolume = createSelectedVolume(secondVolume, firstVolume);

        if (selectedVolume) {
            LINFO("Creating combined volume with dimensions " << selectedVolume->getDimensions() << " using "
                << filteringMode_.get() << " filtering ...");
            copyInputToSelectedVolume(selectedVolume, firstVolume, secondVolume, selectedVolume_.getValue());
        }
    }

    outportFirst_.setData(selectedVolume);
    outportSecond_.setData(selectedVolume);
}

inline bool withinRange(const tgt::vec3& pos, const tgt::vec3& llf, const tgt::vec3& urb) {
    return tgt::hand(tgt::greaterThanEqual(pos, llf)) &&
           tgt::hand(tgt::lessThanEqual(   pos, urb));
}

void VolumeSelectionAndDuplication::copyInputToSelectedVolume(Volume* selectedVolume, const VolumeBase* firstVolume,
                                   const VolumeBase* secondVolume, SelectedVolumeItem itemSelected) const {

    tgtAssert(selectedVolume && firstVolume && secondVolume, "Null pointer passed");

    // compute transformation from voxel coordinates of combined volume
    // to voxel coords of input volumes
    tgt::mat4 combinedToFirst = computeConversionMatrix(selectedVolume, firstVolume);
    tgt::mat4 combinedToSecond = computeConversionMatrix(selectedVolume, secondVolume);

    LDEBUG("Voxel-to-world (First): " << combinedToFirst);
    LDEBUG("Voxel-to-world (Second) " << combinedToSecond);

    // determine which volume is the other (non reference) one
    const VolumeBase* otherVolume = selectedVolume_.isSelected("first") ? secondVolume : firstVolume;

    //
    // voxel-wise combination
    //
    tgt::vec3 dimFirst(firstVolume->getDimensions() - svec3(1));
    tgt::vec3 dimSecond(secondVolume->getDimensions() - svec3(1));
    const bool nearestFiltering = filteringMode_.isSelected("nearest");
    const bool linearFiltering = filteringMode_.isSelected("linear");
    const bool cubicFiltering = filteringMode_.isSelected("cubic");

    const VolumeRAM* v1 = firstVolume->getRepresentation<VolumeRAM>();
    const VolumeRAM* v2 = secondVolume->getRepresentation<VolumeRAM>();
    VolumeRAM* vc = selectedVolume->getWritableRepresentation<VolumeRAM>();

    VRN_FOR_EACH_VOXEL_WITH_PROGRESS(pos, tgt::ivec3(0), selectedVolume->getDimensions(), const_cast<VolumeSelectionAndDuplication*>(this)) {
        for (size_t referenceChannel = 0; referenceChannel < selectedVolume->getNumChannels(); ++referenceChannel) {
            size_t otherChannel = (selectedVolume->getNumChannels() == otherVolume->getNumChannels()) ? referenceChannel : 0;

            // transform sampling pos to coordinate systems of input volumes
            tgt::vec3 posFirst = combinedToFirst*tgt::vec3(pos);
            tgt::vec3 posSecond = combinedToSecond*tgt::vec3(pos);

            // sample input volumes, if transformed voxel position lies inside respective volume
            float valFirst = 0.f;
            float valSecond = 0.f;
            if (nearestFiltering) {
                if (withinRange(posFirst, tgt::vec3::zero, dimFirst))
                    valFirst = v1->getVoxelNormalized(tgt::iround(posFirst), referenceChannel);
                if (withinRange(posSecond, tgt::vec3::zero, dimSecond))
                    valSecond = v2->getVoxelNormalized(tgt::iround(posSecond), otherChannel);
            }
            else if (linearFiltering) {
                if (withinRange(posFirst, tgt::vec3::zero, dimFirst))
                    valFirst = v1->getVoxelNormalizedLinear(posFirst, referenceChannel);
                if (withinRange(posSecond, tgt::vec3::zero, dimSecond))
                    valSecond = v2->getVoxelNormalizedLinear(posSecond, otherChannel);

            }
            else if (cubicFiltering) {
                if (withinRange(posFirst, tgt::vec3::zero, dimFirst))
                    valFirst = v1->getVoxelNormalizedCubic(posFirst, referenceChannel);
                if (withinRange(posSecond, tgt::vec3::zero, dimSecond))
                    valSecond = v2->getVoxelNormalizedCubic(posSecond, otherChannel);
            }
            else {
                LERROR("Unknown filter mode: " << filteringMode_.get());
                return;
            }

            // apply operation to sampled values (note: a switch-block within a volume traversal loop
            // should normally be avoided, however in this case the main operations are way more expensive)
            float result = 0.f;
            switch (itemSelected) {
                case FIRST_VOLUME:
                    result = valFirst;
                    break;
                case SECOND_VOLUME:
                    result = valSecond;
                    break;
                default:
                    LERROR("Unknown operation: " << selectedVolume_.get());
                    return;
            }

            // assign clamped result to combined volume
            vc->setVoxelNormalized(tgt::clamp(result, 0.f, 1.f), pos, referenceChannel);
        }
    } // VRN_FOR_EACH_VOXEL_WITH_PROGRESS
}



void VolumeSelectionAndDuplication::copyInputToSelectedVolumeOnCommonGrid(Volume* selectedVolume, const VolumeBase* firstVolume,
                                               const VolumeBase* secondVolume, SelectedVolumeItem itemSelected) const {
    tgtAssert(selectedVolume && firstVolume && secondVolume, "Null pointer passed");
    tgtAssert(selectedVolume->getDimensions() == firstVolume->getDimensions() &&
              selectedVolume->getDimensions() == secondVolume->getDimensions(), "Volume dimensions mismatch");

    // determine which volume is the other (non reference) one
    const VolumeBase* otherVolume = selectedVolume_.isSelected("first") ? secondVolume : firstVolume;

    const VolumeRAM* v1 = firstVolume->getRepresentation<VolumeRAM>();
    const VolumeRAM* v2 = secondVolume->getRepresentation<VolumeRAM>();
    VolumeRAM* vc = selectedVolume->getWritableRepresentation<VolumeRAM>();

    VRN_FOR_EACH_VOXEL_WITH_PROGRESS(pos, tgt::ivec3(0), selectedVolume->getDimensions(), const_cast<VolumeSelectionAndDuplication*>(this)) {
        for (size_t referenceChannel = 0; referenceChannel < selectedVolume->getNumChannels(); ++referenceChannel) {
            size_t otherChannel = (selectedVolume->getNumChannels() == otherVolume->getNumChannels()) ? referenceChannel : 0;
            float valFirst = v1->getVoxelNormalized(pos, referenceChannel);
            float valSecond = v2->getVoxelNormalized(pos, otherChannel);

            // apply operation to input voxel values
            float result = 0.f;
            switch (itemSelected) {
                case FIRST_VOLUME:
                    result = valFirst;
                    break;
                case SECOND_VOLUME:
                    result = valSecond;
                    break;
                default:
                    LERROR("Unknown operation: " << selectedVolume_.get());
                    return;
            }

            // assign clamped result to combined volume
            //vc->setVoxelNormalized(tgt::clamp(result, 0.f, 1.f), pos, referenceChannel); //FIXME: clamp prevents working with float volumes
            vc->setVoxelNormalized(result, pos, referenceChannel); //FIXME: clamp prevents working with float volumes
        } // for referenceChannel
    } // VRN_FOR_EACH_VOXEL_WITH_PROGRESS
}

Volume* VolumeSelectionAndDuplication::createSelectedVolume(const VolumeBase* refVolume, const VolumeBase* secondVolume) const {
    // compute untransformed bounding box of reference volume
    tgt::vec3 refLLF = refVolume->getLLF();  // all components negative
    tgt::vec3 refURB = refVolume->getURB();  // all components positive
    std::pair<vec3, vec3> refBB = std::pair<vec3, vec3>(refVolume->getLLF(), refVolume->getURB());

    // create bounding box for second volume and transform into reference coordinate system
    MeshGeometry secondBB = MeshGeometry::createCube(secondVolume->getLLF(), secondVolume->getURB());
    tgt::mat4 transformToRef = refVolume->getWorldToPhysicalMatrix();
    transformToRef *= secondVolume->getPhysicalToWorldMatrix();
    secondBB.transform(transformToRef);

    // determine combined bounding box of reference volume and transformed second volume
    tgt::vec3 combinedLLF = refLLF;
    tgt::vec3 combinedURB = refURB;
    for (MeshGeometry::const_iterator face = secondBB.begin(); face != secondBB.end(); ++face) {
        for (FaceGeometry::const_iterator vertex = face->begin(); vertex != face->end(); ++vertex) {
            combinedLLF = tgt::min(combinedLLF, transformToRef * vertex->getCoords());
            combinedURB = tgt::max(combinedURB, transformToRef * vertex->getCoords());
        }
    }
    std::pair<vec3, vec3> combinedBB = std::pair<vec3, vec3>(combinedLLF, combinedURB);
    LDEBUG("Combined BB: " << combinedBB.first << ", " << combinedBB.second);
    tgtAssert(tgt::hand(tgt::lessThanEqual(combinedBB.first, refBB.first)) &&
        tgt::hand(tgt::greaterThanEqual(combinedBB.second, refBB.second)), "Invalid combined bounding box");

    // derive resolution of combined volume from size ratios of combinedBB/refBB and the ref volume's
    // original resolution => spatial resolution of combined volume equals spatial res. of reference volume
    tgt::vec3 scaleFactors = (combinedBB.second - combinedBB.first) / (refBB.second - refBB.first);
    tgtAssert(tgt::hand(tgt::greaterThan(scaleFactors, vec3(0.f))), "invalid scale factors");
    tgt::svec3 combinedDim = svec3(tgt::iceil(vec3(refVolume->getDimensions()) * scaleFactors));
    LDEBUG("Scale factors: " << scaleFactors);
    LDEBUG("Common Dim: " << combinedDim);
    tgtAssert(tgt::hand(tgt::greaterThanEqual(combinedDim, refVolume->getDimensions())),
        "Invalid combined volume dimensions");

    // create combined volume with proper dimensions and spacing
    Volume* selectedVolume = 0;
    try {
        selectedVolume = VolumeOperatorResize::APPLY_OP(refVolume, combinedDim);
    }
    catch (const std::bad_alloc&) {
        LERROR("Failed to create combined volume with dimensions " << combinedDim << " : bad allocation");
        delete selectedVolume;
        return 0;
    }

    // determine correct offset and spacing for combined volume
    selectedVolume->setOffset(combinedBB.first);
    vec3 s = combinedBB.second - combinedBB.first;
    s /= vec3(combinedDim);
    selectedVolume->setSpacing(s);
    selectedVolume->setPhysicalToWorldMatrix(refVolume->getPhysicalToWorldMatrix());

    return selectedVolume;
}

} // namespace
