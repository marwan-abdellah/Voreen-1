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

#include "voreen/qt/widgets/shortcutpreferenceswidget.h"
#include "voreen/qt/widgets/keydetectorwidget.h"
#include "voreen/core/vis/processors/processor.h"
#include "voreen/core/vis/properties/eventproperty.h"
#include "voreen/core/vis/processors/networkevaluator.h"
#include "voreen/qt/widgets/eventpropertywidget.h"

#include <QComboBox>
#include <QLabel>
#include <QVBoxLayout>

namespace voreen {

ShortcutPreferencesWidget::ShortcutPreferencesWidget(NetworkEvaluator* evaluator) 
    : evaluator_(evaluator)
{
    layout_ = new QVBoxLayout(this);
    createWidgets();
}

void ShortcutPreferencesWidget::rebuildWidgets() {
    foreach (EventPropertyWidget* wdt, widgetList_) {
        wdt->hide();
        layout_->removeWidget(wdt);
    }
    createWidgets();
}

void ShortcutPreferencesWidget::createWidgets() {
    std::vector<Processor*> processors = evaluator_->getProcessors();
    std::vector<EventProperty*> props;

    // Get all the Processors which have an EventProperty into a map with the corresponding EventProperties
    for (std::vector<Processor*>::iterator procIter = processors.begin(); procIter != processors.end(); ++procIter) {
        Properties allProps = (*procIter)->getProperties();
        for (Properties::iterator propIter = allProps.begin(); propIter != allProps.end(); ++propIter) {
            EventProperty* eventProp = dynamic_cast<EventProperty*>(*propIter);
            if (eventProp != 0) {
                props.push_back(eventProp);
            }
        }
    }

    for (std::vector<EventProperty*>::iterator iter = props.begin(); iter != props.end(); ++iter) {
        EventPropertyWidget* wdt = new EventPropertyWidget(*iter);
        widgetList_ << wdt;
        layout_->addWidget(wdt);
    }
}


} // namespace