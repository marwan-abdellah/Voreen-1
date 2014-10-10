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

#include "modules/flowreen/flowreenmodule.h"

// processors
#include "processors/flowmagnitudes3d.h"
#include "processors/floworthogonalslicerenderer.h"
#include "processors/flowreenadapter.h"
#include "processors/flowslicerenderer2d.h"
#include "processors/flowslicerenderer3d.h"
#include "processors/flowstreamlinestexture3d.h"
#include "processors/pathlinerenderer3d.h"
#include "processors/streamlinerenderer3d.h"

// I/O
#include "modules/flowreen/io/flowreader.h"

// VolumeOperators
#include "modules/flowreen/datastructures/volumeoperatorintensitymask.h"

namespace voreen {

FlowreenModule::FlowreenModule(const std::string& modulePath)
    : VoreenModule(modulePath)
{
    setID("Flowreen");
    setGuiName("Flowreen");

    // processors
    registerSerializableType(new FlowMagnitudes3D());
    registerSerializableType(new FlowOrthogonalSliceRenderer());
    registerSerializableType(new FlowSliceRenderer2D());
    registerSerializableType(new FlowSliceRenderer3D());
    registerSerializableType(new FlowStreamlinesTexture3D());
    registerSerializableType(new FlowreenAdapter());
    registerSerializableType(new PathlineRenderer3D());
    registerSerializableType(new StreamlineRenderer3D());

    // I/O
    registerVolumeReader(new FlowReader());

    addShaderPath(getModulePath("glsl"));

    INST_SCALAR_TYPES(VolumeOperatorIntensityMask, VolumeOperatorIntensityMaskGeneric)
}

} // namespace
