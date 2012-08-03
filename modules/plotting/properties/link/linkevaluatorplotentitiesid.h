/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2012 University of Muenster, Germany.                        *
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

#ifndef VRN_LINKEVALUATORPLOTENTITIESID_H
#define VRN_LINKEVALUATORPLOTENTITIESID_H

#include "voreen/core/properties/link/linkevaluatorbase.h"
#include "voreen/core/properties/link/linkevaluatoridgeneric.h"

#include "../../datastructures/plotentitysettings.h"

namespace voreen {

class LinkEvaluatorPlotEntitiesId : public LinkEvaluatorIdGeneric<std::vector<PlotEntitySettings> > {
public:
    virtual std::string getClassName() const { return "LinkEvaluatorPlotEntitiesId"; }
    virtual LinkEvaluatorBase* create() const { return new LinkEvaluatorPlotEntitiesId(); }
};

} // namespace

#endif // VRN_LINKEVALUATORPLOTENTITIESID_H