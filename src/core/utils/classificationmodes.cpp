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

#include "voreen/core/utils/classificationmodes.h"
#include "voreen/core/datastructures/volume/volumegl.h"
#include "voreen/core/datastructures/transfunc/transfunc1dkeys.h"
#include "voreen/core/datastructures/transfunc/preintegrationtable.h"
#include "voreen/core/properties/optionproperty.h"

namespace voreen {

using tgt::vec3;

std::string ClassificationModes::getShaderDefineSamplerType(const std::string mode, const TransFunc* tf, const std::string& defineName) {
    if(startsWith(mode, "pre-integrated") && dynamic_cast<const TransFunc1DKeys*>(tf) )
        return "#define " + defineName + " sampler2D\n";
    else if (mode == "transfer-function")
        return tf->getShaderDefines(defineName);
    else
        return "";
}

std::string ClassificationModes::getShaderDefineFunction(const std::string mode, const std::string& defineName) {
    std::string def = "#define " + defineName + "(transFunc, transFuncTex, voxel, lastIntensity) ";

    if (mode == "none")
        def += "vec4(voxel.a);\n";
    else if (mode == "transfer-function")
        def += "applyTF(transFunc, transFuncTex, voxel);\n";
    else if (startsWith(mode, "pre-integrated"))
        def += "applyTFpi(transFunc, transFuncTex, voxel.a, lastIntensity);\n";

    return def;
}

void ClassificationModes::bindTexture(const std::string mode, TransFunc* tf, float samplingStepSize) {
    if (tf) {
        TransFunc1DKeys* tf1d = dynamic_cast<TransFunc1DKeys*>(tf);

        //if(mode == "pre-integrated-fast")
        //    tf1d->getPreIntegrationTable(samplingStepSize, 0, true)->getTexture()->bind();
        /*else*/ if(mode == "pre-integrated-cpu")
            tf1d->getPreIntegrationTable(samplingStepSize, 0, false)->getTexture()->bind();
        else if(mode == "pre-integrated-gpu")
            tf1d->getPreIntegrationTable(samplingStepSize, 0, false, true)->getTexture()->bind();
        else
            tf->bind();
    }
}

bool ClassificationModes::usesTransferFunction(const std::string mode) {
    if(mode == "none")
        return false;
    else if(startsWith(mode, "pre-integrated"))
        return true;
    else if (mode == "transfer-function")
        return true;
    else
        return false;
}

void ClassificationModes::fillProperty(StringOptionProperty* prop) {
    prop->addOption("none", "none");
    prop->addOption("transfer-function", "Transfer Function");
    prop->addOption("pre-integrated-gpu", "Pre-Integrated TF (GPU)");
    prop->addOption("pre-integrated-cpu", "Pre-integrated TF (CPU)");
    //prop->addOption("pre-integrated-fast", "Pre-integrated TF (fast)");
    prop->select("transfer-function");
}

} // namespace
