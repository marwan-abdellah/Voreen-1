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

#include "volumeinversefft.h"

#include "voreen/core/datastructures/volume/volumeatomic.h"
#include "voreen/core/datastructures/volume/operators/volumeoperatorresize.h"
#include "voreen/core/datastructures/geometry/meshlistgeometry.h"

#include <fftw3.h>

using tgt::ivec3;
using tgt::svec3;
using tgt::vec3;

namespace voreen {

const std::string VolumeInvertFFT::loggerCat_("voreen.VolumeInvertFFT");

VolumeInvertFFT::VolumeInvertFFT()
    : CachingVolumeProcessor()
    , inportFirst_(Port::INPORT, "volume.first", "Volume1 Input")
    , inportSecond_(Port::INPORT, "volume.second", "Volume2 Input")
    , outport_(Port::OUTPORT, "outport", "Volume Output", true)
    , enableProcessing_("enabled", "Enable", true)
    , factorC_("factorC", "Factor c", 0.5f, -2.f, 2.f)
{
    addPort(inportFirst_);
    addPort(inportSecond_);
    addPort(outport_);

    factorC_.setTracking(false);

    addProperty(enableProcessing_);
    addProperty(factorC_);
}

VolumeInvertFFT::~VolumeInvertFFT() {}

Processor* VolumeInvertFFT::create() const {
    return new VolumeInvertFFT();
}

void VolumeInvertFFT::process() {
    /*
    tgtAssert(inportFirst_.getData() && inportFirst_.getData()->getRepresentation<VolumeRAM>(), "No input volume");
    tgtAssert(inportSecond_.getData() && inportSecond_.getData()->getRepresentation<VolumeRAM>(), "No input volume");

    // Get the real and the imaginary components of the spectral volume.
    const VolumeBase* firstVolume = inportFirst_.getData();
    const VolumeBase* secondVolume = inportSecond_.getData();

    const int fX = (int) firstVolume->getDimensions().x;
    const int fY = (int) firstVolume->getDimensions().y;
    const int fZ = (int) firstVolume->getDimensions().z;

    const int sX = (int) secondVolume->getDimensions().x;
    const int sY = (int) secondVolume->getDimensions().y;
    const int sZ = (int) secondVolume->getDimensions().z;

    if (!(fX == sX && fY == sY && fZ == sZ))
        return;

    // Volume size
    const int volumeSize = fX * fY * fZ;

    fftwf_complex* spectralVolume = malloc (sizeof(fftwf_complex) * volumeSize);


    const VolumeRAM* vFirstVolume = firstVolume->getRepresentation<VolumeRAM>();
    if (!vFirstVolume)
        return 0;

    const VolumeAtomic<float>* vAtomicFirst = dynamic_cast<const VolumeAtomic<float>*>(vFirstVolume);
    if (!vAtomicFirst)
        return 0;

    VolumeAtomic<float>* output = vAtomicFirst->clone();

    for (int i = 0; i < volumeSize; i++) {
        spectralVolume[i][0] = vAtomicFirst->voxel(i);
        spectralVolume[i][1] = vAtomicFirst->voxel(i);
    }

    Volume* combinedVolume = 0;

    if (!enableProcessing_.get()) {
        outport_.setData(const_cast<VolumeBase*>(inportFirst_.getData()), false);
        return;
    }
    if (firstVolume->getNumChannels() == secondVolume->getNumChannels()) {
        LINFO("Performing channel-wise combination.");
    }
    else {
        if (referenceVolume_.isSelected("first"))
            LINFO("Number of channels of input volumes do not match. Combining each channel of first volume with channel 0 of second volume");
        else
            LINFO("Number of channels of input volumes do not match. Combining each channel of second volume with channel 0 of first volume");
    }


    // optimized combination for volumes that share a common grid in world-space
    if (firstVolume->getDimensions() == secondVolume->getDimensions() &&
        firstVolume->getSpacing() == secondVolume->getSpacing()     &&
        firstVolume->getOffset() == secondVolume->getOffset()     &&
        firstVolume->getPhysicalToWorldMatrix() == secondVolume->getPhysicalToWorldMatrix()) {

        try {
            if (referenceVolume_.isSelected("first"))
                combinedVolume = firstVolume->clone();
            else
                combinedVolume = secondVolume->clone();
        }
        catch (const std::bad_alloc&) {
            LERROR("Failed to create combined volume with dimensions " << firstVolume->getDimensions()
                << " : bad allocation");
        }

        if (combinedVolume) {
            LINFO("Performing optimized combination on common grid with dimensions "
                << firstVolume->getDimensions() << "...");
            combineVolumesOnCommonGrid(combinedVolume, firstVolume, secondVolume, combineFunction_.getValue());
        }
    }
    // standard combination with resampling
    else {
        if (referenceVolume_.isSelected("first"))
            combinedVolume = createCombinedVolume(firstVolume, secondVolume);
        else
            combinedVolume = createCombinedVolume(secondVolume, firstVolume);

        if (combinedVolume) {
            LINFO("Creating combined volume with dimensions " << combinedVolume->getDimensions() << " using "
                << filteringMode_.get() << " filtering ...");
            combineVolumes(combinedVolume, firstVolume, secondVolume, combineFunction_.getValue());
        }
    }

    outport_.setData(combinedVolume);
    */
}

} // namespace
