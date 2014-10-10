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

#include "vectormagnitude.h"
#include "voreen/core/datastructures/volume/volumeram.h"
#include "voreen/core/datastructures/volume/volume.h"
#include "voreen/core/datastructures/volume/operators/volumeoperatormagnitude.h"

namespace voreen {

const std::string VectorMagnitude::loggerCat_("voreen.VectorMagnitude");

VectorMagnitude::VectorMagnitude()
    : CachingVolumeProcessor(),
    inport_(Port::INPORT, "volumehandle.input", "Volume Input"),
    outport_(Port::OUTPORT, "volumehandle.output", "Volume Output", false),
    enableProcessing_("enableProcessing", "Enable")
{
    addPort(inport_);
    addPort(outport_);

    addProperty(enableProcessing_);
}

VectorMagnitude::~VectorMagnitude() {
}

Processor* VectorMagnitude::create() const {
    return new VectorMagnitude();
}

void VectorMagnitude::process() {
    const VolumeBase* inputHandle = inport_.getData();
    const VolumeRAM* inputVolume = inport_.getData()->getRepresentation<VolumeRAM>();
    tgtAssert(inputVolume, "no input volume");

    if (enableProcessing_.get()) {
        Volume* outputHandle = 0;
        VolumeOperatorMagnitude voOpMa;
        tgtAssert(inputVolume->getNumChannels() > 0, "invalid channel count");
        if (inputVolume->getNumChannels() == 3 || inputVolume->getNumChannels() == 4) {
            if (inputVolume->getBitsAllocated() / inputVolume->getNumChannels() > 8)
                outputHandle = voOpMa.apply<uint16_t>(inputHandle);
            else
                outputHandle = voOpMa.apply<uint8_t>(inputHandle);
        }
        else {
            LWARNING("Input volume of format RGB or RGBA expected.");
        }

        // assign computed volume to outport
        if (outputHandle)
            outport_.setData(outputHandle);
        else
            outport_.setData(0);
    }
    else {
        outport_.setData(inputHandle, false);
    }
}

void VectorMagnitude::deinitialize() throw (tgt::Exception) {
    outport_.setData(0);

    VolumeProcessor::deinitialize();
}

}   // namespace
