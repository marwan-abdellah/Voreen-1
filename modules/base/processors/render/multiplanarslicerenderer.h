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

#ifndef VRN_MULTIPLANARSLICERENDERER_H
#define VRN_MULTIPLANARSLICERENDERER_H

#include "voreen/core/processors/volumerenderer.h"
#include "voreen/core/properties/boolproperty.h"
#include "voreen/core/properties/intproperty.h"
#include "voreen/core/properties/cameraproperty.h"
#include "voreen/core/interaction/camerainteractionhandler.h"
#include "voreen/core/datastructures/volume/volumeslicehelper.h"
#include "voreen/core/properties/transfuncproperty.h"

#include "tgt/glmath.h"
#include "tgt/camera.h"

namespace tgt {
    class TextureUnit;
    class Shader;
}

namespace voreen {

class VRN_CORE_API MultiplanarSliceRenderer : public VolumeRenderer {
public:
    MultiplanarSliceRenderer();
    virtual ~MultiplanarSliceRenderer();
    virtual Processor* create() const;

    virtual std::string getClassName() const { return "MultiplanarSliceRenderer"; }
    virtual std::string getCategory() const  { return "Slice Rendering"; }
    virtual CodeState getCodeState() const   { return CODE_STATE_STABLE; }

protected:
    virtual void process();
    virtual void beforeProcess();
    virtual void initialize() throw (tgt::Exception);
    virtual void deinitialize() throw (tgt::Exception);

protected:
    enum TextureMode {
        TEXTURE_2D,
        TEXTURE_3D
    };

    virtual void setDescriptions() {
        setDescription("Renders three orthogonal slices, aligned to the x-, y- and z-axis.");
    }

    virtual void renderSlice(SliceAlignment sliceAlign, int sliceNo, tgt::TextureUnit& texUnit);
    virtual void updatePropertyConfiguration();

    virtual std::string generateHeader();

    /// Recompiles the shader.
    bool rebuildShader();

    RenderPort outport_;
    VolumePort inport_;

    TransFuncProperty transferFunc1_;
    TransFuncProperty transferFunc2_;
    TransFuncProperty transferFunc3_;
    TransFuncProperty transferFunc4_;

    OptionProperty<TextureMode> texMode_;     ///< use 2D slice textures or 3D volume texture?

    BoolProperty renderXYSlice_;
    BoolProperty renderXZSlice_;
    BoolProperty renderYZSlice_;
    IntProperty sliceNumberXY_;
    IntProperty sliceNumberXZ_;
    IntProperty sliceNumberYZ_;
    CameraProperty camProp_;
    CameraInteractionHandler* cameraHandler_;

    tgt::Shader* sliceShader_;

    static const std::string loggerCat_;
};

}   // namespace

#endif
