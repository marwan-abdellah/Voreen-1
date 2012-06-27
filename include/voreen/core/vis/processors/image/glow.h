/**********************************************************************
 *                                                                    *
 * Voreen - The Volume Rendering Engine                               *
 *                                                                    *
 * Copyright (C) 2005-2009 Visualization and Computer Graphics Group, *
 * Department of Computer Science, University of Muenster, Germany.   *
 * <http://viscg.uni-muenster.de>                                     *
 *                                                                    *
 * This file is part of the Voreen software package. Voreen is free   *
 * software: you can redistribute it and/or modify it under the terms *
 * of the GNU General Public License version 2 as published by the    *
 * Free Software Foundation.                                          *
 *                                                                    *
 * Voreen is distributed in the hope that it will be useful,          *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of     *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the       *
 * GNU General Public License for more details.                       *
 *                                                                    *
 * You should have received a copy of the GNU General Public License  *
 * in the file "LICENSE.txt" along with this program.                 *
 * If not, see <http://www.gnu.org/licenses/>.                        *
 *                                                                    *
 * The authors reserve all rights not expressly granted herein. For   *
 * non-commercial academic use see the license exception specified in *
 * the file "LICENSE-academic.txt". To get information about          *
 * commercial licensing please contact the authors.                   *
 *                                                                    *
 **********************************************************************/

#ifndef VRN_GLOW_H
#define VRN_GLOW_H

#include "voreen/core/vis/processors/image/imageprocessor.h"

namespace voreen {

/**
 * Performs a glow.
 */
class Glow : public ImageProcessor {
public:
    /**
     * The Constructor.
     *
     * @param camera The camera from wich we will get information about the current modelview-matrix.
     * @param tc The TextureContainer that will be used to manage TextureUnits for all render-to-texture work done by the PostProcessing.
     */
    Glow();
    virtual const Identifier getClassName() const {return "ImageProcessor.Glow";}
    virtual const std::string getProcessorInfo() const;
    virtual Processor* create() const {return new Glow();}
    void process(LocalPortMapping* portMapping);

    static const Identifier setGlowColor;
    static const Identifier setGlowRange;

protected:

    void renderInternalPass(int source, int pass);
    ColorProp glowColor_;
    IntProp glowRange_;
};


} // namespace voreen

#endif //VRN_GLOW_H