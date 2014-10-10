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

#include <fftw3.h>
#include "volumecreatefftspectrum.h"
#include "volumeformatconversion.h"
#include "voreen/core/datastructures/volume/operators/volumeoperatorconvert.h"
#include "voreen/core/datastructures/volume/operators/volumeoperatorinvert.h"
#include "voreen/core/datastructures/volume/operators/volumeoperatorfft.h"
#include "tgt/vector.h"


namespace voreen {

const std::string VolumeCreateFFTSpectrum::loggerCat_("voreen.base.VolumeCreateFFTSpectrum");

VolumeCreateFFTSpectrum::VolumeCreateFFTSpectrum()
    : CachingVolumeProcessor(),
    inport_(Port::INPORT, "volumehandle.input", "Volume Input"),
    outport_(Port::OUTPORT, "volumehandle.output", "Volume Output", false),
    enableProcessing_("enabled", "Enable", false)
{
    addPort(inport_);
    addPort(outport_);

    addProperty(enableProcessing_);
}

VolumeCreateFFTSpectrum::~VolumeCreateFFTSpectrum() {}

Processor* VolumeCreateFFTSpectrum::create() const {
    return new VolumeCreateFFTSpectrum();
}

void VolumeCreateFFTSpectrum::process() {
    if (!enableProcessing_.get()) {
        LINFO("Spectrum creation disabled");
        // If not enabled by default, stream the input to both the output ports.
        outport_.setData(const_cast<VolumeBase*>(inport_.getData()), false);
    }
    if (enableProcessing_.get() || inport_.hasChanged()) {
        LINFO("Spectrum creation enabled");
        createFFTSpectrum();
    }
}

void VolumeCreateFFTSpectrum::createFFTSpectrum() {

    const VolumeBase* handle = inport_.getData();
    tgtAssert(handle, "Inport has no data");

    if (handle->getRepresentation<VolumeRAM>()) {
        //VolumeRAM* v = handle->getRepresentation<VolumeRAM>()->clone();
        Volume* v =  VolumeOperatorFFT::APPLY_OP(handle);
        outport_.setData(v);
    }
    else {
        outport_.setData(0);
    }
}

}   // namespace
