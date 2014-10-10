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

#include "voreen/core/ports/port.h"

#include "voreen/core/ports/conditions/portcondition.h"
#include "voreen/core/processors/processor.h"
#include "tgt/event/event.h"

#include <sstream>
#include <stdexcept>

namespace voreen {

const std::string Port::loggerCat_("voreen.Port");

Port::Port(PortDirection direction, const std::string& id, const std::string& guiName, bool allowMultipleConnections, Processor::InvalidationLevel invalidationLevel)
    : PropertyOwner(id, guiName)
    , connectedPorts_()
    , processor_(0)
    , direction_(direction)
    , allowMultipleConnections_(allowMultipleConnections)
    , hasChanged_(false)
    , blockEvents_("blockEvents", "Block Events", false)
    , invalidationLevel_(invalidationLevel)
    , isLoopPort_(false)
    , numLoopIterations_(1)
    , currentLoopIteration_(0)
    , initialized_(false)
{
    if (isOutport()) {
        allowMultipleConnections_ = true;

        addProperty(blockEvents_);
        blockEvents_.setGroupID(id);
        setPropertyGroupGuiName(id, (isInport() ? "Inport: " : "Outport: ") + guiName_);
    }
}

Port::~Port() {
    for (size_t i=0; i<conditions_.size(); i++)
        delete conditions_.at(i);
    conditions_.clear();

    if (isInitialized()) {
        LWARNING("~Port() '" << getQualifiedName() << "' has not been deinitialized");
    }

    disconnectAll();
}

Port* Port::create() const {
    return create(static_cast<PortDirection>(0), "", "");
}

void Port::addCondition(PortCondition* condition) {
    tgtAssert(condition, "condition must not be null");
    if (isOutport()) {
        std::string message = "Adding port conditions to an outport is not allowed";
        LERROR(message);
        throw std::invalid_argument(message);
    }

    condition->setCheckedPort(this);
    conditions_.push_back(condition);
}

void Port::setProcessor(Processor* p) {
    tgtAssert(processor_ == 0, "Processor already set!");
    processor_ = p;
}

bool Port::connect(Port* inport) {
    if (testConnectivity(inport)) {
        connectedPorts_.push_back(inport);
        inport->connectedPorts_.push_back(this);
        getProcessor()->invalidate(invalidationLevel_);
        notifyAfterConnectionAdded(inport);
        inport->notifyAfterConnectionAdded(this);
        inport->invalidatePort();
        return true;
    }
    return false;
}

void Port::disconnect(Port* other) {
    tgtAssert(other, "passed null pointer");
    tgtAssert(other != this, "tried to disconnect port from itself");

    for (size_t i = 0; i < connectedPorts_.size(); ++i) {
        if (connectedPorts_[i] == other) {
            notifyBeforeConnectionRemoved(other);
            other->notifyBeforeConnectionRemoved(this);

            connectedPorts_.erase(connectedPorts_.begin() + i);
            other->disconnect(this);
            getProcessor()->invalidate(invalidationLevel_);

            invalidatePort();
            return;
        }
    }
}

void Port::disconnectAll() {
    for (size_t i=0; i < connectedPorts_.size(); ++i) {
        connectedPorts_[i]->disconnect(this);
    }

    connectedPorts_.clear();
}

bool Port::isReady() const {
    if (isInport())
        return isConnected() && checkConditions();
    else
        return isConnected();
}

bool Port::testConnectivity(const Port* inport) const {
    if ((inport == 0) || (inport == this))
        return false;

    // This port must be an outport, the inport must be not an outport and they
    // must not be already connected.
    if ((isInport()) || (inport->isOutport()) || isConnectedTo(inport))
        return false;

    // If the inport must allow multiple connections or may not be already
    // connected to any other port.
    if ((inport->allowMultipleConnections() == false) && inport->isConnected())
        return false;

    // The ports must not belong to the same processor:
    if (getProcessor() == inport->getProcessor())
        return false;

    // The ports must be of the same type:
    if (typeid(*this) != typeid(*inport))
        return false;

    // Only loop ports are allowed to be connected to loop ports
    if (this->isLoopPort() != inport->isLoopPort())
        return false;

    if (detectIllegalLoop(inport))
        return false;

    return true;
}

bool Port::detectIllegalLoop(const Port* inport) const {
    if(isLoopPort())
        return false;

    Processor* curProc = getProcessor();
    Processor* newProc = inport->getProcessor();

    if(curProc == newProc)
        return true;

    std::set<Processor*> regularSuccessors;
    regularSuccessors.insert(curProc);
    regularSuccessors.insert(newProc);

    std::deque<Processor*> q;
    q.push_back(curProc);
    q.push_back(newProc);

    while (!q.empty()) {
        Processor* p = q.front();
        q.pop_front();

        // Get all outports, including co-processor outports
        std::vector<Port*> ports = p->getOutports();
        // TODO I think co-processor "loops" should be legal...? FL
        //std::vector<CoProcessorPort*> coOutports = p->getCoProcessorOutports();
        //for(size_t i = 0; i < coOutports.size(); i++)
            //ports.push_back((Port*)coOutports.at(i));

        for (size_t i = 0; i < ports.size(); ++i) {
            // search only successors which are not part of a loop-port loop
            if(!ports.at(i)->isLoopPort()) {
                std::vector<const Port*> connectedPorts = ports.at(i)->getConnected();
                for (size_t j = 0; j < connectedPorts.size(); ++j) {
                    Processor* connectedProc = connectedPorts.at(j)->getProcessor();
                    bool elemInserted = (regularSuccessors.insert(connectedProc)).second;
                    if(elemInserted)
                        q.push_back(connectedProc);
                    else if(connectedProc == curProc)
                        return true;
                }
            }
        }
    }

    return false;
}

size_t Port::getNumConnections() const {
    return connectedPorts_.size();
}

void Port::setLoopPort(bool isLoopPort) {
    isLoopPort_ = isLoopPort;
}

bool Port::isLoopPort() const {
    return isLoopPort_;
}

int Port::getNumLoopIterations() const {

    if (!isLoopPort()) {
        LWARNING("getNumLoopIterations() called on non loop-port");
        return 0;
    }

    if (isInport())
        return std::max(numLoopIterations_, 1);
    else if (isOutport() && isConnected()) {
        tgtAssert(!getConnected().empty(), "Connected ports vector is empty");
        return getConnected().front()->getNumLoopIterations();
    }
    else {
        return 1;
    }
}

void Port::setNumLoopIterations(int iterations) {

    if (!isLoopPort()) {
        LWARNING("setNumLoopIterations() called on non loop-port");
        return;
    }

    if (numLoopIterations_ != iterations) {
        numLoopIterations_ = iterations;
        if (getProcessor())
            getProcessor()->invalidate(Processor::INVALID_PORTS);
    }
}

int Port::getLoopIteration() const {

    if (!isLoopPort()) {
        LWARNING("getLoopIteration() called on non loop-port");
        return 0;
    }

    if (isInport())
        return currentLoopIteration_;
    else if (isOutport() && isConnected()) {
        tgtAssert(!getConnected().empty(), "Connected ports vector is empty");
        return getConnected().front()->getLoopIteration();
    }
    else {
        return 0;
    }
}

const std::vector<const Port*> Port::getConnected() const {
    std::vector<const Port*> p;
    for(size_t i=0; i<connectedPorts_.size(); i++)
        p.push_back(connectedPorts_[i]);
    return p;
}

bool Port::isConnected() const {
    return (connectedPorts_.size() > 0);
}

bool Port::isConnectedTo(const Port* port) const {
    for (size_t i = 0; i < connectedPorts_.size(); ++i) {
        if (connectedPorts_[i] == port)
            return true;
    }
    return false;
}

void Port::invalidatePort() {
    hasChanged_ = true;
    if (isOutport()) {
        for (size_t i = 0; i <  connectedPorts_.size(); ++i)
             connectedPorts_[i]->invalidatePort();
    }
    else {
        getProcessor()->invalidate(invalidationLevel_);
        forwardData();
    }
}

void Port::invalidate(int inv /*= 1*/) {
    if (getProcessor())
        getProcessor()->invalidate(inv);
}

bool Port::allowMultipleConnections() const {
    return allowMultipleConnections_;
}

Processor* Port::getProcessor() const {
    return processor_;
}

bool Port::isOutport() const {
    return direction_ == OUTPORT;
}

bool Port::isInport() const {
    return direction_ == INPORT;
}

std::string Port::getContentDescription() const {
    std::stringstream strstr;
    strstr << getGuiName() << std::endl
           << "Type: " << getClassName();
    //set empty flag
    if(!hasData())
        strstr << std::endl << "Data: " << "<empty>";
    return strstr.str();
}

std::string Port::getContentDescriptionHTML() const {
    std::stringstream strstr;
    strstr << "<center><font><b>" << getGuiName() << "</b></font></center>"
           << "Type: " << getClassName();
    //set empty flag
    if(!hasData())
        strstr << "<br>" << "Data: " << "&lt;empty&gt;" ;
    return strstr.str();
}

std::string Port::getQualifiedName() const {
    std::string id;
    if (getProcessor())
        id = getProcessor()->getID() + ".";
    id += getID();
    return id;
}

bool Port::hasChanged() const {
    if (isOutport()) {
        LWARNINGC("voreen.port", "Called hasChanged() on outport!");
    }
    return hasChanged_;
}

void Port::setValid() {
    hasChanged_ = false;

    if (isOutport()) {
           LWARNINGC("voreen.port", "Called setValid() on outport!" << getID() );
    }
}

void Port::clear() {
}

bool Port::supportsCaching() const {
    return false;
}

std::string Port::getHash() const {
    return "";
}

void Port::saveData(const std::string& /*path*/) const throw (VoreenException)  {
    throw VoreenException("Port type does not support saving of its data.");
}

void Port::loadData(const std::string& /*path*/) throw (VoreenException) {
    throw VoreenException("Port type does not support loading of its data.");
}

void Port::distributeEvent(tgt::Event* e) {
    if (isOutport()) {
        if (!blockEvents_.get())
            getProcessor()->onPortEvent(e,this);
    } else {
        for (size_t i = 0; i < connectedPorts_.size(); ++i) {
            if (e->isAccepted())
                return;
            connectedPorts_[i]->distributeEvent(e);
        }
    }
}

void Port::toggleInteractionMode(bool interactionMode, void* source) {
    if (isOutport()) {
        for (size_t i = 0; i < connectedPorts_.size(); ++i)
            connectedPorts_[i]->toggleInteractionMode(interactionMode, source);
    } else {
        getProcessor()->toggleInteractionMode(interactionMode, source);
    }
}

void Port::setLoopIteration(int iteration) {
    currentLoopIteration_ = iteration;
    if (currentLoopIteration_ >= getNumLoopIterations())
        LWARNINGC("voreen.Port", "Current loop iteration greater than number of loop iterations");
}

void Port::initialize() throw (tgt::Exception) {

    if (isInitialized()) {
        std::string id;
        if (getProcessor())
            id = getProcessor()->getID() + ".";
        id += getID();
        LWARNING("initialize(): '" << id << "' already initialized");
        return;
    }

    initialized_ = true;
}

void Port::deinitialize() throw (tgt::Exception) {

    if (!isInitialized()) {
        /*std::string id;
        if (getProcessor())
            id = getProcessor()->getName() + ".";
        id += getName();
        LWARNING("deinitialize(): '" << id << "' not initialized"); */
        return;
    }

    initialized_ = false;
}

bool Port::checkConditions() const {
    for (size_t i=0; i<conditions_.size(); i++) {
        if (!conditions_.at(i)->acceptsPortData()) {
            LWARNING("Port condition of '" << getQualifiedName() << "' not met: " << conditions_.at(i)->getDescription());
            return false;
        }
    }
    return true;
}

bool Port::isInitialized() const {
    return initialized_;
}

std::string Port::getDescription() const {
    return description_;
}

void Port::setDescription(std::string desc) {
    description_ = desc;
}

tgt::col3 Port::getColorHint() const {
    return tgt::col3(0, 0, 0);
}

void Port::addForwardPort(Port* port){
    tgtAssert(port->isOutport(), "Only outports can get forwared data!");
    tgtAssert(!port->getClassName().compare(getClassName()),"Forward ports have to be the same type as this class!");

    forwardPorts_.push_back(port);
}

bool Port::removeForwardPort(Port* port){
    for(std::vector<Port*>::iterator it = forwardPorts_.begin(); it != forwardPorts_.end(); it++){
        if(*it == port){
            forwardPorts_.erase(it);
            return true;
        }
    }
    return false;
}

void Port::serialize(XmlSerializer& s) const {
    PropertyOwner::serialize(s);

    s.serialize("direction", direction_);
    s.serialize("portID", id_);
    s.serialize("guiName", guiName_);
    s.serialize("allowMultipleConnections", allowMultipleConnections_);
    s.serialize("invalidationLevel", invalidationLevel_);
}

void Port::deserialize(XmlDeserializer& s) {
    // only deserialize main data fields, if port has been created dynamically via create()
    if (getID().empty()) {
        try {
            int dir;
            s.deserialize("direction", dir);
            direction_ = static_cast<PortDirection>(dir);
            s.deserialize("portID", id_);
            s.deserialize("allowMultipleConnections", allowMultipleConnections_);
            int level;
            s.deserialize("invalidationLevel", level);
            invalidationLevel_ = static_cast<Processor::InvalidationLevel>(level);
        }
        catch (XmlSerializationNoSuchDataException&) {
            s.removeLastError();
        }
    }

    try {
        s.deserialize("guiName", guiName_);
    }
    catch (XmlSerializationNoSuchDataException&) {
        s.removeLastError();
    }

    PropertyOwner::deserialize(s);
}

void Port::notifyAfterConnectionAdded(const Port* connectedPort) {
    std::vector<PortObserver*> observers = Observable<PortObserver>::getObservers();
    for (size_t i=0; i<observers.size(); ++i)
        observers[i]->afterConnectionAdded(this, connectedPort);
}

void Port::notifyBeforeConnectionRemoved(const Port* connectedPort) {
    std::vector<PortObserver*> observers = Observable<PortObserver>::getObservers();
    for (size_t i=0; i<observers.size(); ++i)
        observers[i]->beforeConnectionRemoved(this, connectedPort);
}

} // namespace voreen
