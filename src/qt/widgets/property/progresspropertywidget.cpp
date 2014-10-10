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

#include "voreen/qt/widgets/property/progresspropertywidget.h"

#include "voreen/core/properties/progressproperty.h"
#include "voreen/qt/widgets/customlabel.h"

#include <QProgressBar>
#include <QLabel>
#include <QApplication>
#include <QThread>

namespace voreen {

ProgressPropertyWidget::ProgressPropertyWidget(ProgressProperty* prop, QWidget* parent)
    : QPropertyWidget(prop, parent, false)
    , property_(prop)
    , progressBar_(new QProgressBar())
{

    addWidget(progressBar_);
    progressBar_->setRange(0, 100);
    progressBar_->setValue(0);

    addVisibilityControls();
    QFontInfo fontInfo(font());
}

void ProgressPropertyWidget::updateFromPropertySlot() {
    bool isGuiThread = (QThread::currentThread() == QCoreApplication::instance()->thread());
    if (!isGuiThread)
        return;

    progressBar_->setValue(tgt::iround(property_->getProgress() * 100.f));
    qApp->processEvents();
}


} // namespace
