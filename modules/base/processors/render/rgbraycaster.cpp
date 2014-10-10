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

#include "rgbraycaster.h"

#include "voreen/core/properties/cameraproperty.h"

#include "tgt/textureunit.h"

using tgt::vec3;
using tgt::TextureUnit;

namespace voreen {

RGBRaycaster::RGBRaycaster()
    : VolumeRaycaster()
    , volumePort_(Port::INPORT, "volumehandle.volumehandle", "Volume Input")
    , entryPort_(Port::INPORT, "image.entrypoints", "Entry-points Input", false)
    , exitPort_(Port::INPORT, "image.exitpoints", "Exit-points Input", false)
    , outport_(Port::OUTPORT, "image.output", "Image Output", true)
    , internalRenderPort_(Port::OUTPORT, "internalPort", "Internal Port")
    , shaderProp_("raycast.prg", "Raycasting Shader", "rc_rgb.frag", "passthrough.vert")
    , transferFunc_("transferFunction", "Transfer function", Processor::INVALID_RESULT,
        TransFuncProperty::Editors(TransFuncProperty::INTENSITY))
    , applyColorModulation_("applyColorModulation", "Apply color modulation", false)
    , camera_("camera", "Camera", tgt::Camera(vec3(0.f, 0.f, 3.5f), vec3(0.f, 0.f, 0.f), vec3(0.f, 1.f, 0.f)), true)
{
    addPort(volumePort_);
    addPort(entryPort_);
    addPort(exitPort_);
    addPort(outport_);

    addPrivateRenderPort(internalRenderPort_);

    addProperty(shaderProp_);
    addProperty(transferFunc_);
    addProperty(applyColorModulation_);
    addProperty(camera_);
    camera_.setVisible(false);
}

Processor* RGBRaycaster::create() const {
    return new RGBRaycaster();
}

void RGBRaycaster::initialize() throw (tgt::Exception) {
    VolumeRaycaster::initialize();
    compile();
}

void RGBRaycaster::compile() {
    shaderProp_.setHeader(generateHeader());
    shaderProp_.rebuild();
}

void RGBRaycaster::process() {

    if (!volumePort_.isReady())
        return;

    if (!outport_.isReady())
        return;

    // compile program
    if (getInvalidationLevel() >= Processor::INVALID_PROGRAM)
        compile();
    LGL_ERROR;

    if(volumePort_.hasChanged())
        camera_.adaptInteractionToScene(volumePort_.getData()->getBoundingBox().getBoundingBox());

    // determine render size and activate internal port group
    const bool renderCoarse = interactionMode() && interactionCoarseness_.get() > 1;
    const tgt::svec2 renderSize = (renderCoarse ? (outport_.getSize() / interactionCoarseness_.get()) : outport_.getSize());
    internalRenderPort_.resize(renderSize);
    internalRenderPort_.activateTarget();
    internalRenderPort_.clearTarget();
    LGL_ERROR;

    // initialize shader
    tgt::Shader* raycastPrg = shaderProp_.getShader();
    raycastPrg->activate();

    // set common uniforms used by all shaders
    tgt::Camera cam = camera_.get();
    setGlobalShaderParameters(raycastPrg, &cam, renderSize);
    LGL_ERROR;

    // bind entry/exit param textures
    tgt::TextureUnit entryUnit, entryDepthUnit, exitUnit, exitDepthUnit;
    entryPort_.bindTextures(entryUnit, entryDepthUnit, GL_NEAREST);
    raycastPrg->setUniform("entryPoints_", entryUnit.getUnitNumber());
    raycastPrg->setUniform("entryPointsDepth_", entryDepthUnit.getUnitNumber());
    entryPort_.setTextureParameters(raycastPrg, "entryParameters_");

    exitPort_.bindTextures(exitUnit, exitDepthUnit, GL_NEAREST);
    raycastPrg->setUniform("exitPoints_", exitUnit.getUnitNumber());
    raycastPrg->setUniform("exitPointsDepth_", exitDepthUnit.getUnitNumber());
    exitPort_.setTextureParameters(raycastPrg, "exitParameters_");
    LGL_ERROR;

    // bind volume and pass related information to shader
    std::vector<VolumeStruct> volumeTextures;
    TextureUnit volUnit;
    volumeTextures.push_back(VolumeStruct(
        volumePort_.getData(),
        &volUnit,
        "volume_","volumeStruct_")
    );
    bindVolumes(raycastPrg, volumeTextures, &cam, lightPosition_.get());
    LGL_ERROR;

    // bind transfer function and pass it to the shader
    TextureUnit transferUnit;
    if (transferFunc_.get()) {
        transferUnit.activate();
        transferFunc_.get()->bind();
        transferFunc_.get()->setUniform(raycastPrg, "transferFunc_", "transferFuncTex_", transferUnit.getUnitNumber());
    }

    // perform the actual raycasting by drawing a screen-aligned quad
    renderQuad();

    raycastPrg->deactivate();
    internalRenderPort_.deactivateTarget();
    LGL_ERROR;

    // copy over rendered images from internal buffer to outport,
    // thereby rescaling them to outport dimensions
    rescaleRendering(internalRenderPort_, outport_);
    //outport_.deactivateTarget();

    TextureUnit::setZeroUnit();
    LGL_ERROR;
}

std::string RGBRaycaster::generateHeader() {
    std::string header = VolumeRaycaster::generateHeader();

    header += transferFunc_.get()->getShaderDefines();
    return header;
}


} // namespace voreen
