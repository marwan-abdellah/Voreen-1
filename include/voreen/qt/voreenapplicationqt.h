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

#ifndef VRN_APPLICATIONQT_H
#define VRN_APPLICATIONQT_H

#include "voreen/core/voreenapplication.h"
#include "voreen/core/properties/buttonproperty.h"

//includes for networkeditor
#include "voreen/qt/networkeditor/editor_settings.h"
#include "voreen/core/properties/intproperty.h"
#include "voreen/core/properties/optionproperty.h"
#include "voreen/core/properties/boolproperty.h"
#include "voreen/core/properties/floatproperty.h"

#include "voreen/qt/progressdialog.h"

#include "voreen/qt/voreenqtapi.h"
#include <QString>

class QMainWindow;

namespace voreen {

class VoreenModuleQt;

class VRN_QT_API VoreenApplicationQt : public VoreenApplication {
public:
    VoreenApplicationQt(const std::string& name, const std::string& displayName, const std::string& description,
                        int argc, char** argv, ApplicationFeatures appType = APP_DEFAULT);
    ~VoreenApplicationQt();

    virtual void initialize() throw (VoreenException);
    virtual void deinitialize() throw (VoreenException);

    virtual void initializeGL() throw (VoreenException);
    virtual void deinitializeGL() throw (VoreenException);

    /**
     * Allows access to the global instance of this class.
     */
    static VoreenApplicationQt* qtApp();

    /**
     * Registers a Voreen Qt module.
     */
    void registerQtModule(VoreenModuleQt* module);

    /**
     * Returns all registered Voreen Qt modules.
     */
    const std::vector<VoreenModuleQt*>& getQtModules() const;

    // Returns the VoreenModuleQt specified by the name or 0 if no such module exists
    VoreenModuleQt* getQtModule(const std::string& moduleName) const;

    void setMainWindow(QMainWindow* mainWindow);

    QMainWindow* getMainWindow() const;

    /**
     * Creates a tgt::QtTimer.
     *
     * @param handler The event handler that will be used
     *  for broadcasting the timer events. Must not be null.
     */
    virtual tgt::Timer* createTimer(tgt::EventHandler* handler) const;

    /**
     * Creates a ProgressDialog.
     */
    virtual ProgressDialog* createProgressDialog() const;

    /**
     * Displays a message box.
     */
    virtual void showMessageBox(const std::string& title, const std::string& message, bool error=false) const;

    /**
     * Clears the QSettings used by this application.
     */
    virtual void resetApplicationSettings();

    /**
     * Returns the the scaleProcessorFontSizeProperty value.
     */
    virtual int getProcessorFontScale() const;

    /**
     * Registers a widget to receive touch events.  Default implementation does nothing.
     */
    virtual void sendTouchEventsTo(QWidget *wid) {}

protected:
    virtual void loadModules() throw (VoreenException);

private:
    /**
     * Queries the user whether the Qt application settings should
     * be reset and calls resetApplicationSettings() on positive answer.
     */
    void queryResetApplicationSettings();

    /// Button for resetting the Qt application settings (displayed by the VoreenVE settings dialog).
    ButtonProperty resetApplicationSettingsButton_;

    ///Properties for the network editor
    IntProperty scaleProcessorFontSizeProperty_;
    OptionProperty<NetworkEditorStyles> networkEditorStyleProperty_;
    OptionProperty<NetworkEditorGraphLayouts> networkEditorGraphLayoutsProperty_;
    FloatProperty sugiShiftXProperty_;
    BoolProperty sugiOverlapProperty_;
    BoolProperty sugiMedianProperty_;
    BoolProperty sugiPortFlushProperty_;



    static VoreenApplicationQt* qtApp_;
    QMainWindow* mainWindow_;

    /// indicates that the application settings have to be cleared on application destruction
    bool clearSettings_;

    std::vector<VoreenModuleQt*> qtModules_;

    static const std::string loggerCat_;

};

} // namespace

#endif //VRN_APPLICATIONQT_H
