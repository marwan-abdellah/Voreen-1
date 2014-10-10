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

#include "volumehalfsample.h"
#include "voreen/core/datastructures/volume/volumeram.h"
#include "voreen/core/datastructures/volume/volume.h"
#include "voreen/core/datastructures/volume/operators/volumeoperatorhalfsample.h"

namespace voreen {

const std::string VolumeHalfsample::loggerCat_("voreen.base.VolumeHalfsample");

VolumeHalfsample::VolumeHalfsample()
    : CachingVolumeProcessor(),
    inport_(Port::INPORT, "input", "Volume Input"),
    outport_(Port::OUTPORT, "output", "Volume Output", false),
    enableProcessing_("enableProcessing", "Enable")
{
    addPort(inport_);
    addPort(outport_);

    addProperty(enableProcessing_);
}

VolumeHalfsample::~VolumeHalfsample() {}

Processor* VolumeHalfsample::create() const {
    return new VolumeHalfsample();
}

void VolumeHalfsample::process() {
    const VolumeBase* inputVolume = inport_.getData();

    if (!enableProcessing_.get()) {
        outport_.setData(inputVolume, false);
        return;
    }

    Volume* outputVolume = VolumeOperatorHalfsample::APPLY_OP(inputVolume, this);
    outport_.setData(outputVolume);

}

}   // namespace
