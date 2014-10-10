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

#include "labelmapoverlayimagefilter.h"
#include "voreen/core/datastructures/volume/volumeram.h"
#include "voreen/core/datastructures/volume/volume.h"
#include "voreen/core/datastructures/volume/volumeatomic.h"
#include "voreen/core/ports/conditions/portconditionvolumetype.h"
#include "modules/itk/utils/itkwrapper.h"
#include "voreen/core/datastructures/volume/operators/volumeoperatorconvert.h"
#include "itkImage.h"

#include "itkLabelMapOverlayImageFilter.h"
#include "itkBinaryImageToLabelMapFilter.h"

#include <iostream>

namespace voreen {

const std::string LabelMapOverlayImageFilterITK::loggerCat_("voreen.LabelMapOverlayImageFilterITK");

LabelMapOverlayImageFilterITK::LabelMapOverlayImageFilterITK()
    : ITKProcessor(),
    inport1_(Port::INPORT, "LabelMap"),
    inport2_(Port::INPORT, "FeatureImage"),
    outport1_(Port::OUTPORT, "OutputImage"),
    numObjects_("numObjects", "Number of Objects", 0, 0, 1000),
    opacity_("opacity", "Opacity", 0.5f, 0.f, 1.f)
{
    addPort(inport1_);
    PortConditionLogicalOr* orCondition1 = new PortConditionLogicalOr();
    orCondition1->addLinkedCondition(new PortConditionVolumeTypeUInt8());
    orCondition1->addLinkedCondition(new PortConditionVolumeTypeInt8());
    orCondition1->addLinkedCondition(new PortConditionVolumeTypeUInt16());
    orCondition1->addLinkedCondition(new PortConditionVolumeTypeInt16());
    orCondition1->addLinkedCondition(new PortConditionVolumeTypeUInt32());
    orCondition1->addLinkedCondition(new PortConditionVolumeTypeInt32());
    inport1_.addCondition(orCondition1);
    addPort(inport2_);
    PortConditionLogicalOr* orCondition2 = new PortConditionLogicalOr();
    orCondition2->addLinkedCondition(new PortConditionVolumeTypeUInt8());
    orCondition2->addLinkedCondition(new PortConditionVolumeTypeInt8());
    orCondition2->addLinkedCondition(new PortConditionVolumeTypeUInt16());
    orCondition2->addLinkedCondition(new PortConditionVolumeTypeInt16());
    orCondition2->addLinkedCondition(new PortConditionVolumeTypeUInt32());
    orCondition2->addLinkedCondition(new PortConditionVolumeTypeInt32());
    inport2_.addCondition(orCondition2);
    addPort(outport1_);

    addProperty(numObjects_);
    numObjects_.setWidgetsEnabled(false);
    addProperty(opacity_);

}

Processor* LabelMapOverlayImageFilterITK::create() const {
    return new LabelMapOverlayImageFilterITK();
}

template<class T, class S>
void LabelMapOverlayImageFilterITK::labelMapOverlayImageFilterITK() {


    typedef itk::Image<T, 3> InputImageType1;
    typedef itk::Image<S, 3> InputImageType2;
    typedef itk::Image<itk::CovariantVector<T,3>, 3> OutputImageType1;

    typename InputImageType1::Pointer p1 = voreenToITK<T>(inport1_.getData());
    typename InputImageType2::Pointer p2 = voreenToITK<S>(inport2_.getData());

    typedef itk::BinaryImageToLabelMapFilter<InputImageType1> LabelMapType;

    //Alternative:
    //typedef itk::LabelImageToLabelMapFilter<InputImageType1> LabelMapType;

    typename LabelMapType::Pointer labelMap = LabelMapType::New();
    labelMap->SetInput(p1);
    labelMap->Update();
    numObjects_.set(labelMap->GetOutput()->GetNumberOfLabelObjects());

    //Filter define
    typedef typename itk::LabelMapOverlayImageFilter<typename LabelMapType::OutputImageType, InputImageType2, OutputImageType1> FilterType;
    typename FilterType::Pointer filter = FilterType::New();

    filter->SetInput(labelMap->GetOutput());
    filter->SetFeatureImage(p2);

    filter->SetOpacity(opacity_.get());


    observe(filter.GetPointer());

    try
    {
        filter->Update();

    }
    catch (itk::ExceptionObject &e)
    {
        LERROR(e);
    }


    Volume* outputVolume1 = 0;
    outputVolume1 = ITKVec3ToVoreenVec3Copy<T>(filter->GetOutput());

    if (outputVolume1) {
        outport1_.setData(outputVolume1);
        transferTransformation(inport1_.getData(), outputVolume1);
    }
    else
        outport1_.setData(0);



}



void LabelMapOverlayImageFilterITK::process() {
    const VolumeBase* inputHandle1 = inport1_.getData();
    const VolumeRAM* inputVolume1 = inputHandle1->getRepresentation<VolumeRAM>();

    if (dynamic_cast<const VolumeRAM_UInt8*>(inputVolume1))  {
        volumeTypeSwitch1<uint8_t>();
    }
    else if (dynamic_cast<const VolumeRAM_Int8*>(inputVolume1))  {
        volumeTypeSwitch1<int8_t>();
    }
    else if (dynamic_cast<const VolumeRAM_UInt16*>(inputVolume1))  {
        volumeTypeSwitch1<uint16_t>();
    }
    else if (dynamic_cast<const VolumeRAM_Int16*>(inputVolume1))  {
        volumeTypeSwitch1<int16_t>();
    }
    else if (dynamic_cast<const VolumeRAM_UInt32*>(inputVolume1))  {
        volumeTypeSwitch1<uint32_t>();
    }
    else if (dynamic_cast<const VolumeRAM_Int32*>(inputVolume1))  {
        volumeTypeSwitch1<int32_t>();
    }
    else {
        LERROR("Inputformat of Volume 1 is not supported!");
    }

}

template <class T>
void LabelMapOverlayImageFilterITK::volumeTypeSwitch1() {
    const VolumeBase* inputHandle2 = inport2_.getData();
    const VolumeRAM* inputVolume2 = inputHandle2->getRepresentation<VolumeRAM>();

    if (dynamic_cast<const VolumeRAM_UInt8*>(inputVolume2))  {
        labelMapOverlayImageFilterITK<T, uint8_t>();
    }
    else if (dynamic_cast<const VolumeRAM_Int8*>(inputVolume2))  {
        labelMapOverlayImageFilterITK<T, int8_t>();
    }
    else if (dynamic_cast<const VolumeRAM_UInt16*>(inputVolume2))  {
        labelMapOverlayImageFilterITK<T, uint16_t>();
    }
    else if (dynamic_cast<const VolumeRAM_Int16*>(inputVolume2))  {
        labelMapOverlayImageFilterITK<T, int16_t>();
    }
    else if (dynamic_cast<const VolumeRAM_UInt32*>(inputVolume2))  {
        labelMapOverlayImageFilterITK<T, uint32_t>();
    }
    else if (dynamic_cast<const VolumeRAM_Int32*>(inputVolume2))  {
        labelMapOverlayImageFilterITK<T, int32_t>();
    }
    else {
        LERROR("Inputformat of Volume 2 is not supported!");
    }

}


}   // namespace
