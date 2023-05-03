import os

import SimpleITK as sitk
from functionsRegistration import resliceVolumeITK
from functionsRegistration import resample_img

#------------------INPUTS-------------------

workFolder = "Registration/"

group = "E11.0_1"
Downsampled = False
scale = 0.7
sampleObjectiveGroup= "Feb12_E115_2"

RegTypeRigid        = "LM38RigidReg_itk"
RegTypeAvgAffine    = "MeanNNTransform"

#PIL TIFFs with dimension 1, but MATLAB generate TIFFs with 0.35...
#Set according to the file. Landmarks are always using voxel size = 1
constDivide = 1

sampleNames = ["Feb3_E11_4_v2","Feb12_E115_1",\
               "Feb12_E115_2","Feb12_E115_8",\
               "May10_E11_2"]


for sampleName in sampleNames:

    print("---------------------------------------")
    print(sampleName)

    #--------------OBJECTIVE-------------------------
    pathVolumeObjective = workFolder + group + "/" + sampleObjectiveGroup + "_Objective" + "/" \
        + sampleObjectiveGroup + "_Tissues_Objective.tiff"

    #-------------- Moving ------------------------------

    folder =        workFolder+ group +"/"+sampleName+"/"
    regTypeFolderAffine = workFolder + group + "/" + RegTypeRigid + "/"

    try:
        os.mkdir(regTypeFolderAffine)
    except OSError as error:
        print(error)

    folderTransformation = regTypeFolderAffine + sampleName + "_to_Objective/"
    folderOutput = regTypeFolderAffine + sampleName + "_to_Objective_OtherVolumes/"

    try:
        os.mkdir(folderTransformation)
    except OSError as error:
        print(error)

    pathRigidTransformationTXT = folderTransformation + sampleName + "_transf_to_Objective.txt"
    pathRigidTransformationTFM = folderTransformation + sampleName + "_transf_to_Objective.tfm"

    pathRigidTransformationInvTXT = folderTransformation + sampleName + "_transf_inverse_to_Objective.txt"
    pathRigidTransformationInvTFM = folderTransformation + sampleName + "_transf_inverse_to_Objective.tfm"

    try:
        os.mkdir(folderOutput)
    except OSError as error:
        print(error)

    regTypeFolderAffine = workFolder + group + "/" + RegTypeRigid + "/" + RegTypeAvgAffine + "/"

    if Downsampled:
        regTypeFolderAffine = workFolder + group + "/" + RegTypeRigid + "/" "Downsampled" + "/"  + RegTypeAvgAffine + "/"

    print(regTypeFolderAffine)
    folderAverageAffine = regTypeFolderAffine + sampleName + "/"
    print(folderAverageAffine)
    pathAvgAffineTransformationTXT = folderAverageAffine + sampleName + "_Tissues_Masked_to_MeanTransform.txt"
    print(pathAvgAffineTransformationTXT)
    folderOutputAffine = regTypeFolderAffine + sampleName + "_OtherVolumes/"

    try:
        os.mkdir(folderOutputAffine)
    except OSError as error:
        print(error)

    #-------------- Moving ------------------------------

    pathVolumeFixedCorrected = folderTransformation + "FixedCorrected.tiff"

    pathMovingVolumes = [ \
                        #folder + sampleName + "_Cells_Masked.tiff", \
                        #folder + sampleName + "_ProliferatingCells_Masked.tiff", \
                        #folder + sampleName + "_Mesen.tiff", \
                        folder + sampleName + "_TissueSlices_192_processed_3000_Masked.tiff", \
                        ]

    pathMovedVolumes = [\
                        #folderOutput + sampleName + "_Cells_Masked_Reg.tiff", \
                        #folderOutput + sampleName + "_ProliferatingCells_Masked_Reg.tiff", \
                        #folderOutput + sampleName + "_MesenReg.tiff", \
                        folderOutput + sampleName + "_TissueSlices_192_processed_3000_Masked_Reg.tiff", \
                        ]

    pathMovedAvgAffineVolumes = [ \
                        #folderOutputAffine + sampleName + "_Cells_Masked_AvgAffine.tiff", \
                        #folderOutputAffine + sampleName + "_ProliferatingCells_Masked_AvgAffine.tiff", \
                        #folderOutputAffine + sampleName + "_MesenReg_AvgAffine.tiff", \
                        folderOutputAffine + sampleName + "_TissueSlices_192_processed_3000_Masked_AvgAffine.tiff", \
                        ]

    #Get spacing from volumes
    #Read volumes

    moving_im = sitk.ReadImage(pathMovingVolumes[0])
    spacingEmbryoMoving = moving_im.GetSpacing()
    #print("spacingEmbryoMoving" + str(spacingEmbryoMoving))

    fixed_im = sitk.ReadImage(pathVolumeObjective)
    print("Read ITK successful", pathVolumeObjective)

    extent2 = fixed_im
    spacingEmbryoFixed = fixed_im.GetSpacing()
    #print("spacingEmbryoFixed" + str(spacingEmbryoFixed))
    #Registration of landmarks to the fixed

    landmarkTransformITK = sitk.ReadTransform(pathRigidTransformationTXT)
    #landmarkTransformITKInverse = sitk.ReadTransform(pathTransformationInvTXT)

    nVolumes = len(pathMovingVolumes)

    for i in range(nVolumes):

        pathMovingVolume = pathMovingVolumes[i]
        pathResultVolume = pathMovedVolumes[i]
        pathResultAvgAffineVolume = pathMovedAvgAffineVolumes[i]

        print("------------------RIGID REGISTRATION BASED ON LMs------------------------------")

        resliceVolumeITK(pathMovingVolume, pathResultVolume, pathVolumeObjective, pathVolumeFixedCorrected,
                         landmarkTransformITK, constDivide = constDivide)

        print("------------------AVERAGE AFFINE REGISTRATION BASED ON TISSUES------------------------------")

        moving_image = sitk.ReadImage(pathResultVolume, sitk.sitkFloat32)
        moving_image.SetSpacing([constDivide, constDivide, constDivide])

        imageSize = moving_image.GetSize()
        YInic = int(imageSize[1] * 0.35)
        XInic = int(imageSize[0] * 0.15)
        ZInic = int(imageSize[2] * 0.15)
        if group == "E9.5":
            # For E9.5
            YInic = int(imageSize[1] * 0.6)
            XInic = int(imageSize[0] * 0.25)
            ZInic = int(imageSize[2] * 0.25)

        if group == "E11.5":
            # For E11.5
            YInic = 0  # int(imageSize[1] * 0)
            XInic = 0  # int(imageSize[0] * 0)
            ZInic = 0  # int(imageSize[2] * 0)

        moving_image_cropped = moving_image[XInic:(imageSize[0] - XInic), YInic:imageSize[1], ZInic:(imageSize[2] - ZInic)]
        #moving_image_cropped_Size = moving_image_cropped.GetSize()

        if Downsampled:
            original_spacing = moving_image_cropped.GetSpacing()
            out_spacing = [x / scale for x in original_spacing]
            moving_image_cropped = resample_img(moving_image_cropped, out_spacing, is_label=True)
            output_size = moving_image_cropped.GetSize()
            #I have to set the spacing to constDivide
            moving_image_cropped.SetSpacing([constDivide, constDivide, constDivide])
            print('scale downsample: '+str(scale))
            print('output size: '+str(output_size[0]))

        parameterMap = sitk.ReadParameterFile(pathAvgAffineTransformationTXT)
        #print(parameterMap)
        #sitk.PrintParameterMap(parameterMap)
        #resultImageAffine = sitk.Transformix(moving_image_cropped, parameterMap)

        transformixImageFilter = sitk.TransformixImageFilter()
        transformixImageFilter.SetTransformParameterMap(parameterMap)
        transformixImageFilter.SetMovingImage(moving_image_cropped)
        transformixImageFilter.LogToConsoleOff()
        resultImageAffine = transformixImageFilter.Execute()

        #sitk.WriteImage(sitk.Cast(moving_image_cropped, sitk.sitkUInt8), pathResultAvgAffineVolume+'_HowItWasCropped.tiff')
        sitk.WriteImage(sitk.Cast(resultImageAffine, sitk.sitkUInt8), pathResultAvgAffineVolume)
        print("Avg affine transformation done!")



