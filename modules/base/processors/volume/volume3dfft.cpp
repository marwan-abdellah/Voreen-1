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

#include "volume3dfft.h"

#include <fftw3.h>

#include "voreen/core/datastructures/volume/volumeatomic.h"
#include "voreen/core/datastructures/volume/operators/volumeoperatorresize.h"
#include "voreen/core/datastructures/geometry/meshlistgeometry.h"

using tgt::ivec3;
using tgt::svec3;
using tgt::vec3;

namespace voreen {

const std::string Volume3DFFT::loggerCat_("voreen.Volume3DFFT");

Volume3DFFT::Volume3DFFT()
    : CachingVolumeProcessor()
    , inport_(Port::INPORT, "inport", "Volume Input")
    , outportFirst_(Port::OUTPORT, "outport.first", "Volume 1 Output", true)
    , outportSecond_(Port::OUTPORT, "outport.second", "Volume 2 Output", true)
    , enableProcessing_("enabled", "Enable", true)
{
    addPort(inport_);
    addPort(outportFirst_);
    addPort(outportSecond_);

    enableProcessing_.onChange(CallMemberAction<Volume3DFFT>(this, &Volume3DFFT::forceUpdate));
    addProperty(enableProcessing_);
}

Volume3DFFT::~Volume3DFFT() {}

Processor* Volume3DFFT::create() const {
    return new Volume3DFFT();
}

void Volume3DFFT::process() {
    // Getting the input volumes.
    tgtAssert(inport_.getData() && inport_.getData()->getRepresentation<VolumeRAM>(), "No input volume");

    const VolumeBase* inputVolume = inport_.getData();

    // Allocate the output volume and copy the input volume to it.
    Volume* outputVolume1 = 0;
    Volume* outputVolume2 = 0;

    outputVolume1 = inputVolume->clone();
    outputVolume2 = inputVolume->clone();

    // If not processingm just by-pass the first volume.
    if (!enableProcessing_.get() && forceUpdate_) {
        outportFirst_.setData(const_cast<VolumeBase*>(inport_.getData()));
        outportSecond_.setData(const_cast<VolumeBase*>(inport_.getData()));

        return;
    } else {
        // Copy the input volume to the outport ports.
        copyInputVolumeToOutputVolume(outputVolume1, inputVolume);
        copyInputVolumeToOutputVolume(outputVolume2, inputVolume);

        computeForwardFFT(outputVolume1, outputVolume2, inputVolume);
    }

    outportFirst_.setData(outputVolume1);
    outportSecond_.setData(outputVolume2);
}

inline bool withinRange(const tgt::vec3& pos, const tgt::vec3& llf, const tgt::vec3& urb) {
    return tgt::hand(tgt::greaterThanEqual(pos, llf)) &&
           tgt::hand(tgt::lessThanEqual(   pos, urb));
}

void Volume3DFFT::copyInputVolumeToOutputVolume(Volume* output,
                                                      const VolumeBase* input) {

    tgtAssert(input && output, "Null pointer passed");
    tgtAssert(input->getDimensions() == output->getDimensions(), "Volume dimensions mismatch");

    const VolumeRAM* volumeInput = input->getRepresentation<VolumeRAM>();
    VolumeRAM* volumeOutput = output->getWritableRepresentation<VolumeRAM>();

    VRN_FOR_EACH_VOXEL_WITH_PROGRESS(pos, tgt::ivec3(0), input->getDimensions(), const_cast<Volume3DFFT*>(this)) {
        for (size_t referenceChannel = 0; referenceChannel < input->getNumChannels(); ++referenceChannel) {
            float value = volumeInput->getVoxelNormalized(pos, referenceChannel);

            // apply operation to input voxel values
            float result = value;

            // assign clamped result to combined volume
            volumeOutput->setVoxelNormalized(result, pos, referenceChannel); //FIXME: clamp prevents working with float volumes
        } // for referenceChannel
    } // VRN_FOR_EACH_VOXEL_WITH_PROGRESS
}


void Volume3DFFT::computeForwardFFT(Volume* realOutput, Volume* imaginaryOutput, const VolumeBase* input) {

    tgtAssert(input && realOutput && imaginaryOutput, "Null pointer passed");
    tgtAssert(input->getDimensions() == output->getDimensions(), "Volume dimensions mismatch");

    const VolumeRAM* volumeInput = input->getRepresentation<VolumeRAM>();
    VolumeRAM* voreal = realOutput->getWritableRepresentation<VolumeRAM>();
    VolumeRAM* voimaginary = imaginaryOutput->getWritableRepresentation<VolumeRAM>();

    // Allocate the spectral volume
    const uint sx = input->getDimensions().x;
    const uint sy = input->getDimensions().y;
    const uint sz = input->getDimensions().z;
    const uint sVolume = sx * sy * sz;
    fftwf_complex* spectrum = new fftwf_complex[sVolume];

    // Fill in the spectral volume
    VRN_FOR_EACH_VOXEL_WITH_PROGRESS(pos, tgt::ivec3(0), input->getDimensions(), const_cast<Volume3DFFT*>(this)) {
        for (size_t referenceChannel = 0; referenceChannel < input->getNumChannels(); ++referenceChannel) {
            const int index = pos.x + sx * (pos.y + sy * pos.z);
            spectrum[index][0] = volumeInput->getVoxelNormalized(pos, referenceChannel);
            spectrum[index][1] = 0.f;

        } // for referenceChannel
    } // VRN_FOR_EACH_VOXEL_WITH_PROGRESS

    // Compute the 3d fft
    fftwf_plan fftPlan = fftwf_plan_dft_3d(sx, sy, sz, spectrum, spectrum,
                   FFTW_FORWARD, FFTW_ESTIMATE);
    fftwf_execute(fftPlan);

    // Fillout the output volumes.
    VRN_FOR_EACH_VOXEL_WITH_PROGRESS(pos, tgt::ivec3(0), input->getDimensions(), const_cast<Volume3DFFT*>(this)) {
        for (size_t referenceChannel = 0; referenceChannel < input->getNumChannels(); ++referenceChannel) {
            const int index = pos.x + sx * (pos.y + sy * pos.z);
            voreal->setVoxelNormalized(spectrum[index][0], pos, referenceChannel); //FIXME: clamp prevents working with float volumes
            voimaginary->setVoxelNormalized(spectrum[index][1], pos, referenceChannel); //FIXME: clamp prevents working with float volumes

        } // for referenceChannel
    } // VRN_FOR_EACH_VOXEL_WITH_PROGRESS

    // Destroy the fft plan.
    fftwf_destroy_plan(fftPlan);

    // Free the spectral array
    delete [] spectrum;
}

void Volume3DFFT::fillNullVolume(Volume* output) {

    tgtAssert(output, "Null pointer passed");

    VolumeRAM* volumeOutput = output->getWritableRepresentation<VolumeRAM>();

    VRN_FOR_EACH_VOXEL_WITH_PROGRESS(pos, tgt::ivec3(0), output->getDimensions(), const_cast<Volume3DFFT*>(this)) {
        for (size_t referenceChannel = 0; referenceChannel < output->getNumChannels(); ++referenceChannel) {

            // apply operation to input voxel values
            float result = 0.0;

            // assign clamped result to combined volume
            volumeOutput->setVoxelNormalized(result, pos, referenceChannel); //FIXME: clamp prevents working with float volumes
        } // for referenceChannel
    } // VRN_FOR_EACH_VOXEL_WITH_PROGRESS
}

void Volume3DFFT::forceUpdate() {
    forceUpdate_ = true;
}

} // namespace
