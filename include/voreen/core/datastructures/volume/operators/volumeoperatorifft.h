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

#ifndef VRN_VOLUMEOPERATORIFFT_H
#define VRN_VOLUMEOPERATORIFFT_H

#include "voreen/core/datastructures/volume/volumeoperator.h"

#include <fftw3.h>

namespace voreen {

// Base class, defines interface for the operator (-> apply):
class VRN_CORE_API VolumeOperatorIFFTBase : public UnaryVolumeOperatorBase {
public:
    virtual Volume* apply(const VolumeBase* volume) const = 0;
};

// Generic implementation:
template<typename T>
class VolumeOperatorIFFTGeneric : public VolumeOperatorIFFTBase {
public:
    virtual Volume* apply(const VolumeBase* volume) const;
    //Implement isCompatible using a handy macro:
    IS_COMPATIBLE
};

template<typename T>
Volume* VolumeOperatorIFFTGeneric<T>::apply(const VolumeBase* vh) const {
    const VolumeRAM* v = vh->getRepresentation<VolumeRAM>();
    if (!v)
        return 0;

    const VolumeAtomic<T>* va = dynamic_cast<const VolumeAtomic<T>*>(v);
    if (!va)
        return 0;

    VolumeAtomic<T>* output = va->clone();

    // Get volume dimensions.
    tgt::svec3 volDim = vh->getDimensions();
    const uint volumeWidth = volDim.x;
    const uint volumeHeight = volDim.y;
    const uint volumeDepth = volDim.z;
    const uint volumeSize = volumeWidth * volumeHeight * volumeDepth;

    // Allocate the complex (spectral) volume.
    fftwf_complex* spectralVolume = (fftwf_complex*)
            malloc (sizeof(fftwf_complex) * volumeSize);

    // Fill the spectral array
    for (size_t i = 0; i < volumeSize; i++) {
        spectralVolume[i][0] = float(va->voxel(i));
        spectralVolume[i][1] = 0.f;
    }

    fftwf_plan fftPlan = fftwf_plan_dft_3d(volumeWidth, volumeHeight, volumeDepth,
                                      spectralVolume, spectralVolume,
                                      FFTW_BACKWARD, FFTW_ESTIMATE);
    fftwf_execute(fftPlan);

    for (size_t i = 0; i < volumeSize; i++) {
        output->voxel(i) = spectralVolume[i][0];
    }

//    VRN_FOR_EACH_VOXEL(index, tgt::svec3(0, 0, 0), output->getDimensions()) {
//        output->voxel(index) = va->voxel(index) / 2;
//    }

    return new Volume(output, vh);
}

typedef UniversalUnaryVolumeOperatorGeneric<VolumeOperatorIFFTBase> VolumeOperatorIFFT;

} // namespace

#endif // VRN_VOLUMEOPERATORIFFT_H
