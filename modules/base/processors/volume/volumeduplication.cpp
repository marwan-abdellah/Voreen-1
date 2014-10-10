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

#include "volumeduplication.h"

#include "voreen/core/datastructures/volume/volumeatomic.h"
#include "voreen/core/datastructures/volume/operators/volumeoperatorresize.h"
#include "voreen/core/datastructures/geometry/meshlistgeometry.h"

using tgt::ivec3;
using tgt::svec3;
using tgt::vec3;

namespace voreen {

const std::string VolumeDuplication::loggerCat_("voreen.VolumeDuplication");

VolumeDuplication::VolumeDuplication()
    : CachingVolumeProcessor()
    , inport_(Port::INPORT, "inport", "Volume Input")
    , outportFirst_(Port::OUTPORT, "outport.first", "Volume 1 Output", true)
    , outportSecond_(Port::OUTPORT, "outport.second", "Volume 2 Output", true)
    , enableProcessing_("enabled", "Enable", true)
{
    addPort(inport_);
    addPort(outportFirst_);
    addPort(outportSecond_);

    enableProcessing_.onChange(CallMemberAction<VolumeDuplication>(this, &VolumeDuplication::forceUpdate));
    addProperty(enableProcessing_);
}

VolumeDuplication::~VolumeDuplication() {}

Processor* VolumeDuplication::create() const {
    return new VolumeDuplication();
}

void VolumeDuplication::process() {
    // Getting the input volumes.
    tgtAssert(inportFirst_.getData() && inportFirst_.getData()->getRepresentation<VolumeRAM>(), "No input volume");
    tgtAssert(inportSecond_.getData() && inportSecond_.getData()->getRepresentation<VolumeRAM>(), "No input volume");

    const VolumeBase* inputVolume = inport_.getData();

    // Allocate the output volume and copy the input volume to it.
    Volume* outputVolume1 = 0;
    Volume* outputVolume2 = 0;

    outputVolume1 = inputVolume->clone();
    outputVolume2 = inputVolume->clone();

    // If not processingm just by-pass the first volume.
    if (!enableProcessing_.get() && forceUpdate_) {
        // Zero
        outportFirst_.setData(0);
        outportSecond_.setData(0);

        return;
    } else {
        // Copy the input volume to the outport ports.
        copyInputVolumeToOutputVolume(outputVolume1, inputVolume);
        copyInputVolumeToOutputVolume(outputVolume2, inputVolume);
    }

    outportFirst_.setData(outputVolume1);
    outportSecond_.setData(outputVolume2);
}

inline bool withinRange(const tgt::vec3& pos, const tgt::vec3& llf, const tgt::vec3& urb) {
    return tgt::hand(tgt::greaterThanEqual(pos, llf)) &&
           tgt::hand(tgt::lessThanEqual(   pos, urb));
}

void VolumeDuplication::copyInputVolumeToOutputVolume(Volume* output,
                                                      const VolumeBase* input) {

    tgtAssert(input && output, "Null pointer passed");
    tgtAssert(input->getDimensions() == output->getDimensions(), "Volume dimensions mismatch");

    const VolumeRAM* volumeInput = input->getRepresentation<VolumeRAM>();
    VolumeRAM* volumeOutput = output->getWritableRepresentation<VolumeRAM>();

    VRN_FOR_EACH_VOXEL_WITH_PROGRESS(pos, tgt::ivec3(0), input->getDimensions(), const_cast<VolumeDuplication*>(this)) {
        for (size_t referenceChannel = 0; referenceChannel < input->getNumChannels(); ++referenceChannel) {
            float value = volumeInput->getVoxelNormalized(pos, referenceChannel);

            // apply operation to input voxel values
            float result = value;

            // assign clamped result to combined volume
            volumeOutput->setVoxelNormalized(result, pos, referenceChannel); //FIXME: clamp prevents working with float volumes
        } // for referenceChannel
    } // VRN_FOR_EACH_VOXEL_WITH_PROGRESS
}

void VolumeDuplication::fillNullVolume(Volume* output) {

    tgtAssert(output, "Null pointer passed");

    VolumeRAM* volumeOutput = output->getWritableRepresentation<VolumeRAM>();

    VRN_FOR_EACH_VOXEL_WITH_PROGRESS(pos, tgt::ivec3(0), output->getDimensions(), const_cast<VolumeDuplication*>(this)) {
        for (size_t referenceChannel = 0; referenceChannel < output->getNumChannels(); ++referenceChannel) {

            // apply operation to input voxel values
            float result = 0.0;

            // assign clamped result to combined volume
            volumeOutput->setVoxelNormalized(result, pos, referenceChannel); //FIXME: clamp prevents working with float volumes
        } // for referenceChannel
    } // VRN_FOR_EACH_VOXEL_WITH_PROGRESS
}

void VolumeDuplication::forceUpdate() {
    forceUpdate_ = true;
}

} // namespace
