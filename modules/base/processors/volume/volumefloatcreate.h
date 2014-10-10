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

#ifndef VRN_VOLUMEFLOATCREATE_H
#define VRN_VOLUMEFLOATCREATE_H

#include "volumecreatebase.h"

#include "voreen/core/properties/optionproperty.h"
#include "voreen/core/properties/intproperty.h"
#include "voreen/core/properties/buttonproperty.h"
#include "voreen/core/properties/boolproperty.h"

#include <string>

namespace voreen {

class VRN_CORE_API VolumeFloatCreate : public VolumeCreateBase {
public:
    VolumeFloatCreate();
    Processor* create() const { return new VolumeFloatCreate; }

    std::string getClassName() const  { return "VolumeFloatCreate";      }
    std::string getCategory() const   { return "Volume Processing"; }
    CodeState getCodeState() const    { return CODE_STATE_STABLE;   }
    bool usesExpensiveComputation() const { return true; }

protected:
    virtual void setDescriptions() {
        setDescription("Generates an 32-bit dataset with cubic dimensions.");
    }

    void process();

private:
    void createCornell(const tgt::ivec3& dimensions, VolumeRAM_Float *target);
    void createCube(const tgt::ivec3& dimensions, VolumeRAM_Float* target);
    void createBlobs(const tgt::ivec3& dimensions, VolumeRAM_Float* target);
    void createBlobs2(const tgt::ivec3& dimensions, VolumeRAM_Float* target);
    void createBlobs3(const tgt::ivec3& dimensions, VolumeRAM_Float* target);
    void createSphere(const tgt::ivec3& dimensions, VolumeRAM_Float* target);
    void createDoubleSphere(const tgt::ivec3& dimensions, VolumeRAM_Float* target);
    void createTorus(const tgt::ivec3& dimensions, VolumeRAM_Float* target);
    void createDoubleTorus(const tgt::ivec3& dimensions, VolumeRAM_Float* target);
    void createDoublePartialTorus(const tgt::ivec3& dimensions, VolumeRAM_Float* target);
    void createBumpySphere(const tgt::ivec3& dimensions, VolumeRAM_Float* target);
    void createSphereCoord(const tgt::ivec3& dimensions, VolumeRAM_4xUInt8* target);
    void createSynth(const tgt::ivec3& dimensions, VolumeRAM_Float* target);
    void createCloud(const tgt::ivec3& dimensions, VolumeRAM_Float* target);
    void createAOTestBox(const tgt::ivec3& dimensions, VolumeRAM_Float* target);
    void createShadowTestVolume(const tgt::ivec3& dimensions, VolumeRAM_Float* target);
    void createAorticArch(const tgt::ivec3& dimensions, VolumeRAM_Float* target);
    void createRandomShapes(const tgt::ivec3& dimensions, VolumeRAM_Float* target);
    void createDECT(const tgt::ivec3& dimensions, VolumeRAM_Float* target);
    void createStack(const tgt::ivec3& dimensions, VolumeRAM_Float* target);

    VolumePort outport_;
    unsigned int currentSeed_;

    static const std::string loggerCat_; ///< category used in logging

    StringOptionProperty operation_;
    IntProperty dimension_;
    ButtonProperty regenerate_;
    IntProperty numShapes_;
    BoolProperty keepCurrentShapes_;

    IntProperty numSubdivisions_;
};

}   //namespace

#endif
