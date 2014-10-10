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

#ifndef VRN_SLICEVIEWER_H
#define VRN_SLICEVIEWER_H

#include "voreen/core/processors/volumerenderer.h"
#include "voreen/core/properties/eventproperty.h"
#include "voreen/core/properties/optionproperty.h"
#include "voreen/core/properties/boolproperty.h"
#include "voreen/core/properties/intproperty.h"
#include "voreen/core/properties/matrixproperty.h"
#include "voreen/core/interaction/mwheelnumpropinteractionhandler.h"
#include "voreen/core/datastructures/volume/volumeslicehelper.h"
#include "voreen/core/properties/transfuncproperty.h"
#include "voreen/core/datastructures/volume/volumeslicehelper.h"

namespace voreen {

/**
 * Performs slice rendering of a single or multiple slices
 * along one of the three main axis of the volume.
 */
class VRN_CORE_API SliceViewer : public VolumeRenderer {

public:
    SliceViewer();
    virtual ~SliceViewer();
    virtual Processor* create() const;

    virtual std::string getClassName() const    { return "SliceViewer";      }
    virtual std::string getCategory() const     { return "Slice Rendering";  }
    virtual CodeState getCodeState() const      { return CODE_STATE_STABLE;  }

protected:

    virtual void initialize() throw (tgt::Exception);
    virtual void deinitialize() throw (tgt::Exception);

    virtual void beforeProcess();
    virtual void process();
    virtual void afterProcess();

    virtual void interactionModeToggled();

    virtual void adjustPropertiesToInput();

    /// Generates the header for the shader depending on the choice of features to be used.
    virtual std::string generateHeader();

    /// Recompiles the shader.
    bool rebuildShader();

    /**
     * Adapts the min/max ranges of the respective properties to the
     * dimensions of the currently connected volume and adjusts the property visibility.
     */
    void updatePropertyConfiguration();

    /**
     * Adapts the volume coord permutation and the states
     * of the slice properties when the slice alignment has changed.
     */
    void onSliceAlignmentChange();

    /**
     * Renders a GL_LINE_LOOP with white color and the current slice's
     * dimensions in order to visualize its boundary.
     */
    void renderSliceBoundaries() const;

    /**
     * Displays the cursor position in volume coordinates, the data value under the cursor,
     * and/or the current slice number. What data is actually shown, is determined by the
     * respective properties.
     *
     * @note requires freetype (module 'fontrendering')
     */
    void renderInfoTexts() const;

    /**
     * Draws the scale legend in the lower left corner.
     */
    void renderLegend();

    /**
     * Tiny helper function returning true when numSlicesPerRow_ and
     * numSlicesPerCol_ are set to one.
     */
    bool singleSliceMode() const;

    /**
     * Converts the given screen coordinates into the corresponding voxel
     * coordinates for the slice.
     */
    tgt::vec3 screenToVoxelPos(tgt::ivec2 screenPos) const;

    /**
     * Generates a matrix that maps from viewport coordinates
     * to volume coordinates.
     *
     * @note This matrix is only valid in single slice mode,
     *  since for multi slice mode, the mapping from viewport
     *  to volume coordinates cannot be expressed by a linear function.
     */
    tgt::mat4 generatePickingMatrix() const;

    /**
     * Is called by the respective event properties on mouse move or
     * mouse click events and saves the current cursor position.
     */
    void mouseLocalization(tgt::MouseEvent* e);

    /**
     * Is called by the mouseEventShift_ property when a slice shift
     * event has occurred.
     */
    void shiftEvent(tgt::MouseEvent* e);

protected:
    enum TextureMode {
        TEXTURE_2D,
        TEXTURE_3D
    };

    virtual void setDescriptions() {
        setDescription("Displays 2D slices along one of the three main axis of the volume. Multiple slices can be viewed simultaneously.");
    }

    VolumePort inport_;
    RenderPort outport_;

    TransFuncProperty transferFunc1_;
    TransFuncProperty transferFunc2_;
    TransFuncProperty transferFunc3_;
    TransFuncProperty transferFunc4_;

    /// Property containing the available alignments: xy (axial), xz (coronal), yz (sagittal)
    OptionProperty<SliceAlignment> sliceAlignment_;
    IntProperty sliceIndex_;                ///< Property containing the currently selected slice
    IntProperty numGridRows_;               ///< Property containing the current row count of the displayed grid
    IntProperty numGridCols_;               ///< Property containing the current column count of the displayed grid
    BoolProperty selectCenterSliceOnInputChange_;   ///< if true, the center slice index is selected when the input volume changes

    IntProperty mouseXCoord_;               ///< Property containing the current x slice position of the mouse cursor
    IntProperty mouseYCoord_;               ///< Property containing the current y slice position of the mouse cursor
    IntProperty mouseZCoord_;               ///< Property containing the current z slice position of the mouse cursor

    BoolProperty renderSliceBoundaries_;    ///< Determines whether to render the slice boundaries
    FloatVec4Property boundaryColor_;       ///< Color to be used for rendering of the slice boundary
    IntProperty boundaryWidth_;           ///< line width of the boundary

    StringOptionProperty showCursorInfos_;  ///< Determines whether information about the cursor position is to be shown
    BoolProperty showSliceNumber_;          ///< Determines whether the slice number is to displayed on each slice
    IntProperty fontSize_;                  ///< Font size to be used for info texts
    BoolProperty showScaleLegend_;             ///< true, if the distance legend should be rendered

    FloatVec2Property voxelOffset_;         ///< The 2D voxel position at which the view area is shifted
    FloatProperty zoomFactor_;              ///< Specifies the current slice zoom factor: the standard value of 1.0
                                            ///  causes the whole slice to be displayed

    FloatMat4Property pickingMatrix_;       ///< Read-only (generated) matrix mapping from viewport to
                                            ///  volume coordinates.

    OptionProperty<TextureMode> texMode_;   ///< use 2D slice textures or 3D volume texture?
    IntProperty sliceLevelOfDetail_;        ///< level of detail used for 2D slice extraction
    IntProperty interactionLevelOfDetail_;  ///< level of detail during user interaction
    IntProperty sliceExtractionTimeLimit_;  ///< timelimit in milliseconds for 2D slice extraction
    IntProperty sliceCacheSize_;            ///< size of the 2D slice cache (0 means no cache).

    EventProperty<SliceViewer>* mouseEventPress_;
    EventProperty<SliceViewer>* mouseEventMove_;
    EventProperty<SliceViewer>* mouseEventShift_;

    MWheelNumPropInteractionHandler<int> mwheelCycleHandler_;
    MWheelNumPropInteractionHandler<float> mwheelZoomHandler_;

    tgt::Shader* sliceShader_;

    VolumeSliceCache sliceCache_;       ///< Cache for slices created in 2D texture mode.

    tgt::ivec3 voxelPosPermutation_;    ///< permutation of voxel coordinates according to current alignment
    tgt::vec2 sliceLowerLeft_;          ///< lower left of the whole slice plate in viewport coordinates
    tgt::vec2 sliceSize_;               ///< size of a single slice in viewport coordinates
    tgt::mat4 textureMatrix_;           ///< matrix that is currently applied to texture coordinates
    bool sliceComplete_;                ///< is set in process() and specifies whether the current slice has been created in full LOD,
                                        ///< if not, an invalidation is triggered in afterProcess().

    static const std::string fontName_; ///< path and name of the font used for text-rendering

    static const std::string loggerCat_;

private:
    mutable tgt::ivec2 mousePosition_;          ///< Current mouse position
    mutable bool mouseIsPressed_;               ///< Is a mouse button currently pressed?
    mutable tgt::ivec3 lastPickingPosition_;    ///< Mouse position during previous interaction

};

} // namespace

#endif // VRN_SLICEVIEWER_H
