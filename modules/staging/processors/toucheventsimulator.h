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

#ifndef VRN_TOUCHEVENTSIMULATOR_H
#define VRN_TOUCHEVENTSIMULATOR_H

#include "voreen/core/processors/processor.h"
#include "voreen/core/ports/renderport.h"
#include "voreen/core/processors/imageprocessor.h"
#include "tgt/textureunit.h"
#include "voreen/core/properties/eventproperty.h"
#include "tgt/event/touchevent.h"

#include <deque>

namespace voreen {

/**
 * Simulates Touch Events from mouse events when pressing modifier STRG.
 */
class TouchEventSimulator : public ImageProcessor {

public:

    TouchEventSimulator();

    ~TouchEventSimulator();

    virtual Processor* create() const;

    virtual std::string getClassName() const;

    virtual std::string getCategory() const;

protected:

    virtual void setDescriptions();

    /**
     * Bypass image.
     */
    virtual void process();

    /**
     * Catches mouse events when modifier STRG is pressed and produces corresponding touch events.
     */
    void onMouseEvent(tgt::MouseEvent* e);

private:

    RenderPort inport_;           ///< inport receives an image to be passed through
    RenderPort outport_;          ///< outport used to pass through the received image

    EventProperty<TouchEventSimulator>* mouseEventProp_; ///< property for catching the mouse events
};

} // namespace

#endif // VRN_TOUCHEVENTSIMULAROR_H
