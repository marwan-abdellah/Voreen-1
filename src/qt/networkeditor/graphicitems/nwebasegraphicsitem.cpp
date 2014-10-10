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

#include "voreen/qt/networkeditor/graphicitems/nwebasegraphicsitem.h"
#include "voreen/qt/networkeditor/graphicitems/tooltips/tooltipbasegraphicsitem.h"
#include "voreen/qt/networkeditor/styles/nwestyle_base.h"

#include <QGraphicsSceneHoverEvent>
#include <QAction>
#include <QMenu>

namespace voreen {

NWEBaseGraphicsItem::NWEBaseGraphicsItem(NetworkEditor* nwe)
    : QGraphicsObject(0)
    , networkEditor_(nwe)
    , toolTipItem_(0)
    , paintHasBeenInitialized_(false)
    , contextMenuActionList_()
{
    tgtAssert(nwe,"NWEBaseGraphicsItem: No NetworkEditor passed!");
}

NWEBaseGraphicsItem::~NWEBaseGraphicsItem() {
}

//---------------------------------------------------------------------------------------------------------------
//                  getter and setter
//---------------------------------------------------------------------------------------------------------------
NetworkEditor* NWEBaseGraphicsItem::getNetworkEditor() const{
    return networkEditor_;
}

NetworkEditorLayer NWEBaseGraphicsItem::currentLayer() const {
    return networkEditor_->getCurrentLayer();
}

NetworkEditorCursorMode NWEBaseGraphicsItem::currentCursorMode() const {
    return networkEditor_->getCurrentCursorMode();
}

NWEStyle_Base* NWEBaseGraphicsItem::currentStyle() const{
    return networkEditor_->getCurrentStyle();
}

bool NWEBaseGraphicsItem::currentToolTipMode() const{
    return networkEditor_->getCurrentToolTipMode();
}

void NWEBaseGraphicsItem::prepareGeometryChange(){
    QGraphicsItem::prepareGeometryChange();
}

void NWEBaseGraphicsItem::setToolTipGraphicsItem(ToolTipBaseGraphicsItem* tooltip) {
    toolTipItem_ = tooltip;
}

ToolTipBaseGraphicsItem* NWEBaseGraphicsItem::getToolTipGraphicsItem() const {
    return toolTipItem_;
}

bool NWEBaseGraphicsItem::isPaintInitialized() const {
    return paintHasBeenInitialized_;
}

void NWEBaseGraphicsItem::resetPaintInitialization() {
    paintHasBeenInitialized_ = false;
}

//---------------------------------------------------------------------------------------------------------------
//                  paint
//---------------------------------------------------------------------------------------------------------------
void NWEBaseGraphicsItem::paint(QPainter* painter, const QStyleOptionGraphicsItem* option, QWidget* widget){
    if(!paintHasBeenInitialized_) {
        initializePaintSettings();
        layoutChildItems();
        paintHasBeenInitialized_ = true;
    }
    NWEItemSettings setting = beforePaint(painter,option);
    mainPaint(painter,option,widget,setting);
}

NWEItemSettings NWEBaseGraphicsItem::beforePaint(QPainter* painter, const QStyleOptionGraphicsItem* option){
    return currentStyle()->getNWEItemSettings(this,option);
}

//---------------------------------------------------------------------------------------------------------------
//                  context menu
//---------------------------------------------------------------------------------------------------------------
void NWEBaseGraphicsItem::addActionsToContextMenu(QMenu* menu) {
    if(acceptHoverEvents())
        menu->addActions(contextMenuActionList_);
}

void NWEBaseGraphicsItem::setContextMenuActions() {
    foreach(QAction* a, contextMenuActionList_)
        delete a;
    contextMenuActionList_.clear();
}
//---------------------------------------------------------------------------------------------------------------
//                  events
//---------------------------------------------------------------------------------------------------------------
void NWEBaseGraphicsItem::hoverEnterEvent (QGraphicsSceneHoverEvent* event){
    if(toolTipItem_ && currentToolTipMode()){
        QPointF p = event->scenePos();
        toolTipItem_->setToolTipTimerTriggertMousePosition(p);
        toolTipItem_->startTimer();
    }
    QGraphicsItem::hoverEnterEvent(event);
}

void NWEBaseGraphicsItem::mouseMoveEvent (QGraphicsSceneMouseEvent* event){
    if(!currentToolTipMode() && toolTipItem_) {
        if(toolTipItem_->isVisible()){
            toolTipItem_->setVisible(false);
        } else {
            toolTipItem_->stopTimer();
        }
    }
    QGraphicsItem::mouseMoveEvent(event);
}

void NWEBaseGraphicsItem::hoverLeaveEvent (QGraphicsSceneHoverEvent* event){
    if(toolTipItem_){
        toolTipItem_->stopTimer();
        toolTipItem_->setVisible(false);
    }
    QGraphicsItem::hoverLeaveEvent(event);
}


} // namespace voreen
