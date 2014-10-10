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

#include "deprecatedmodule.h"

#include "processors/buttonoverlayprocessor.h"
#include "processors/geometryclippingwidget.h"
#include "processors/textseriessource.h"
#include "processors/volumeseriessource.h"
#include "processors/rawtexturesource.h"
#include "processors/rawtexturesave.h"
#include "processors/targettotexture.h"
#include "processors/texturetotarget.h"
#include "processors/volumenormalization.h"
#include "processors/canny.h"

#include "io/philipsusvolumereader.h"

#include "operators/volumeoperatornormalize.h"

namespace voreen {

DeprecatedModule::DeprecatedModule(const std::string& modulePath)
    : VoreenModule(modulePath)
{
    setID("Deprecated");
    setGuiName("Deprecated");

    addShaderPath(getModulePath("glsl"));

    registerSerializableType(new ButtonOverlayProcessor());
    registerSerializableType(new GeometryClippingWidget());
    registerSerializableType(new TextSeriesSource());
    registerSerializableType(new VolumeSeriesSource());
    registerSerializableType(new RawTextureSource());
    registerSerializableType(new RawTextureSave());
    registerSerializableType(new TargetToTexture);
    registerSerializableType(new TextureToTarget());
    registerSerializableType(new VolumeNormalization());
    registerSerializableType(new Canny());

    registerVolumeReader(new PhilipsUSVolumeReader());

    INST_SCALAR_TYPES(VolumeOperatorNormalize, VolumeOperatorNormalizeGeneric)
    //INST_VECTOR_TYPES(VolumeOperatorNormalize, VolumeOperatorNormalizeGeneric)
}

} // namespace
