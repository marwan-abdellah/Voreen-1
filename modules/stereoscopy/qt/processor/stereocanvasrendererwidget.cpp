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

#include "stereocanvasrendererwidget.h"

#include "voreen/core/network/networkevaluator.h"
#include "voreen/core/utils/voreenpainter.h"
#include "voreen/qt/voreenapplicationqt.h"

#include <QGridLayout>
#include <QGLFormat>

namespace voreen {

StereoCanvasRendererWidget::StereoCanvasRendererWidget(QWidget* parent, StereoCanvasRenderer* canvasRenderer)
    : CanvasRendererWidget(parent, canvasRenderer)
{
}

StereoCanvasRendererWidget::~StereoCanvasRendererWidget() {
}

void StereoCanvasRendererWidget::initialize() throw (VoreenException){
    QProcessorWidget::initialize();

    StereoCanvasRenderer* canvasRenderer = dynamic_cast<StereoCanvasRenderer*>(processor_);
    tgtAssert(canvasRenderer, "StereoCanvasRenderer expected");

    if (!VoreenApplication::app()) {
        LERRORC("voreen.qt.StereoCanvasRendererWidget", "VoreenApplication not instantiated");
        throw VoreenException("VoreenApplication not instantiated");
    }
    NetworkEvaluator* evaluator = VoreenApplication::app()->getNetworkEvaluator(processor_);
    if (!evaluator) {
        LERRORC("voreen.qt.StereoCanvasRendererWidget", "No evaluator assigned to VoreenApplication");
        throw VoreenException("No evaluator assigned to VoreenApplication");
    }

    //check, ich quad-buffering is supported
    QGLFormat stereoTestFormat;
    stereoTestFormat.setStereo(true);
    if(stereoTestFormat.stereo()) {
        canvasWidget_ = new tgt::QtCanvas("", tgt::ivec2(getSize().x, getSize().y), tgt::GLCanvas::RGBADDQ, this, true, 0);
    } else {
        canvasWidget_ = new tgt::QtCanvas("", tgt::ivec2(getSize().x, getSize().y), tgt::GLCanvas::RGBADD, this, true, 0);
        OptionProperty<StereoCanvasRenderer::StereoMode>* prop = dynamic_cast<OptionProperty<StereoCanvasRenderer::StereoMode>* >(canvasRenderer->getProperty("stereoModeProp"));
        tgtAssert(prop, "No stereoModeProp found");
        if(!prop->removeOption("quadbuffer")) {
            tgtAssert(false, "Option \"quadbuffer\" does not exist!");
            LERRORC("voreen.qt.StereoCanvasRendererWidget", "Option \"quadbuffer\" does not exist!");
        }
    }
    VoreenApplicationQt::qtApp()->sendTouchEventsTo(canvasWidget_);
    canvasWidget_->setMinimumSize(64, 64);
    canvasWidget_->setMouseTracking(true); // for receiving mouse move events without a pressed button

    QGridLayout* layout = new QGridLayout;
    layout->setContentsMargins(0, 0, 0, 0);
    layout->addWidget(canvasWidget_, 0, 0);
    setLayout(layout);

    //show();

    VoreenPainter* painter = new VoreenPainter(canvasWidget_, evaluator, canvasRenderer);
    canvasWidget_->setPainter(painter);
    painter->initialize();

    canvasRenderer->setCanvas(canvasWidget_);
    initialized_ = true;
}

} //namespace voreen
