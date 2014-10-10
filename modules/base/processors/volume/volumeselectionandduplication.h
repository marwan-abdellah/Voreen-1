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

#ifndef VRN_VOLUMESELECTIONANDDUPLICATION_H
#define VRN_VOLUMESELECTIONANDDUPLICATION_H

#include "voreen/core/processors/volumeprocessor.h"
#include "voreen/core/properties/boolproperty.h"
#include "voreen/core/properties/floatproperty.h"
#include "voreen/core/properties/optionproperty.h"

namespace voreen {

class VRN_CORE_API VolumeSelectionAndDuplication : public CachingVolumeProcessor {
public:
    VolumeSelectionAndDuplication();
    ~VolumeSelectionAndDuplication();
    virtual Processor* create() const;

    virtual std::string getClassName() const      { return "VolumeSelectionAndDuplication";     }
    virtual std::string getCategory() const       { return "Volume Processing"; }
    virtual CodeState getCodeState() const        { return CODE_STATE_STABLE;   }
    virtual bool usesExpensiveComputation() const { return true; }

protected:
    virtual void setDescriptions() {
        setDescription("Combines two volumes based on a selectable function.");
    }

    virtual void process();

private:
    /// Voxel-wise selection item.
    enum SelectedVolumeItem {
        FIRST_VOLUME,
        SECOND_VOLUME
    };
    friend class OptionProperty<SelectedVolumeItem>;


    /**
     * Copies the input volumes and writes the result to combinedVolume (which is assumed to be already created),
     * by transforming the input volume's coordinates systems to the coordinates system of the combined volume.
     */
    void copyInputToSelectedVolume(Volume* combinedVolume, const VolumeBase* firstVolume,
        const VolumeBase* secondVolume, SelectedVolumeItem item) const;

    /**
     * Copies the input volumes and writes the result to combinedVolume (which is assumed to be already created),
     * without coordinate transformation. This is much faster than the transformation-based combination,
     * but only possible for volumes that share a common grid in world space.
     */
    void copyInputToSelectedVolumeOnCommonGrid(Volume* combinedVolume, const VolumeBase* firstVolume,
        const VolumeBase* secondVolume, SelectedVolumeItem item) const;


    /// Creates a selected (empty) volume from the two input volumes in world space,
    /// or 0 in case the combined volume could not be created due to bad allocation.
    Volume* createSelectedVolume(const VolumeBase* refVolume, const VolumeBase* secondVolume) const;

    VolumePort inportFirst_;
    VolumePort inportSecond_;
    VolumePort outportFirst_;
    VolumePort outportSecond_;

    BoolProperty enableProcessing_;
    OptionProperty<SelectedVolumeItem> selectedVolume_;
    StringOptionProperty filteringMode_;

    static const std::string loggerCat_; ///< category used in logging
};


} // namespace

#endif // VRN_VOLUMESELECTIONANDDUPLICATION_H
