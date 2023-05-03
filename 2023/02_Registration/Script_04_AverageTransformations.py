import os
import SimpleITK as sitk
from functionsLandmarks import transformAndSaveLMsElastix
from functionsRegistration import resample_img
import copy
#from ElastixParameterFile import ElastixParameterFile

workFolder = "Registration/"
group = "E11.0_1"
constDivide = 1
Downsampled = False
sampleObjectiveGroup= "Feb12_E115_2"


RegType = "MeanNNTransform"
previousReg = "LM38RigidReg_itk"
NNReg = "AffineGroupwiseNN"

sampleNames = ["Feb3_E11_4_v2","Feb3_E11_4_v2_flipped","Feb12_E115_1","Feb12_E115_1_flipped",\
               "Feb12_E115_2","Feb12_E115_2_flipped","Feb12_E115_8","Feb12_E115_8_flipped",\
               "May10_E11_2","May10_E11_2_flipped"]


for sampleName in sampleNames:
    print("----------Moving-----------------------------")
    print(sampleName)

    folderPreviousReg   = workFolder + group + "/" + previousReg + "/"
    folderInput         = folderPreviousReg + sampleName + "_to_Objective/"
    regTypeFolder       = workFolder +  group + "/" + previousReg + "/" + RegType + "/"

    if Downsampled:
        folderPreviousReg = workFolder + group + "/" + previousReg + "/" "Downsampled" + "/"
        folderInput = folderPreviousReg + sampleName + "_to_Objective/"
        regTypeFolder = workFolder + group + "/" + previousReg + "/" "Downsampled" + "/" + RegType + "/"

    try:
        os.mkdir(regTypeFolder)
    except OSError as error:
        print(error)

    folderOutput = regTypeFolder + sampleName + "/"

    try:
        os.mkdir(folderOutput)
    except OSError as error:
        print(error)

    pathLMsFaceMoving   = folderInput + sampleName + "_Tissues_Masked_face_to_Objective.pts"
    pathLMsNeckMoving   = folderInput + sampleName + "_Tissues_Masked_neck_to_Objective.pts"
    pathLMsGMMoving     = folderInput + sampleName + "_Tissues_Masked_38lms_to_Objective.pts"

    if Downsampled:
        pathLMsFaceMoving = folderInput + sampleName + "_Tissues_Masked_face_to_Objective_downsampled.pts"
        pathLMsNeckMoving = folderInput + sampleName + "_Tissues_Masked_neck_to_Objective_downsampled.pts"
        pathLMsGMMoving = folderInput + sampleName + "_Tissues_Masked_38lms_to_Objective_downsampled.pts"

    #pathLMsFaceCorrected    = folderOutput + sampleName + "_Tissues_Masked_face_Corrected.pts"
    #pathLMsNeckCorrected    = folderOutput + sampleName + "_Tissues_Masked_neck_Corrected.pts"
    pathLMsGMModifiedForTransform      = folderOutput + sampleName + "_Tissues_Masked_38lms_ModifiedForTransform.pts"
    pathLMsGMCorrected = folderOutput + sampleName + "_Tissues_Masked_38lms_Corrected.pts"

    #e9.5
    #pathLMsFaceModifiedForTransform = folderOutput + sampleName + "_Tissues_Masked_face_ModifiedForTransform.pts"
    #pathLMsFaceCorrected = folderOutput + sampleName + "_Tissues_Masked_face_Corrected.pts"

    pathLMsFaceMoved = folderOutput + sampleName + "_face_to_MeanTransform.pts"
    pathLMsGMMoved = folderOutput + sampleName + "_38lms_to_MeanTransform.pts"

    TransformationDirTXT    = folderOutput  + sampleName + "_Tissues_Masked_to_MeanTransform.txt"
    TransformationDirHDF5   = folderOutput  + sampleName + "_Tissues_Masked_to_MeanTransform.hdf5"
    TransformationDirTFM    = folderOutput  + sampleName + "_Tissues_Masked_to_MeanTransform.tfm"

    infoFile                = folderOutput + sampleName + "_info.txt"

    MovingImage = folderInput + sampleName + "_Tissues_Reg.tiff"

    if Downsampled:
        MovingImage = folderInput + sampleName + "_downsampled_Tissues_Reg.tiff"

    moving_image = sitk.ReadImage(MovingImage, sitk.sitkFloat32)

    moving_image.SetSpacing([constDivide, constDivide, constDivide])
    print(moving_image.GetSpacing())

    imageSize = moving_image.GetSize()

    # Between E10.0 and E11.0
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
        YInic = 0 #int(imageSize[1] * 0)
        XInic = 0 #int(imageSize[0] * 0)
        ZInic = 0 #int(imageSize[2] * 0)

    moving_image_cropped = moving_image[XInic:(imageSize[0]-XInic), YInic:imageSize[1], ZInic:(imageSize[2]-ZInic)]

    moving_image_cropped_Size = moving_image_cropped.GetSize()

    del moving_image

    conversionFiles = []
    TransformationDir = []

    #Write information file
    strInfo = "Original Dimensions\n"
    strInfo = strInfo + "X: " + str(imageSize[0]) + "\nY: " + str(imageSize[1]) + "\nZ: " + str(imageSize[2]) + "\n\n"
    strInfo = strInfo + "Crop:\n"
    strInfo = strInfo + "XInic: " + str(XInic) + "\nYInic: " + str(YInic) + "\nZInic: " + str(ZInic) + "\n\n"
    strInfo = strInfo + "Final Dimensions:\n"
    strInfo = strInfo + "X: " + str(moving_image_cropped_Size[0]) + "\nY: " + str(moving_image_cropped_Size[1]) + "\nZ: " + str(moving_image_cropped_Size[2]) + "\n"
    f = open(infoFile, "w")
    f.write(strInfo)
    f.close()

    for sampleObjectiveGroup in sampleNames:

        #if not (sampleName == sampleObjectiveGroup):
        if (1 == 1):

            # sampleName =    "Dec2_E10_16_flipped"
            print("----------Img to Img---------------")
            print(sampleObjectiveGroup)

            folder = workFolder + group + "/" + previousReg + "/" + NNReg + "/"

            if Downsampled:
                folder = workFolder + group + "/" + previousReg + "/" "Downsampled" + "/"  + NNReg + "/"

            folderNN = folder + sampleName + "_to_" + sampleObjectiveGroup
            TransformationDir = folderNN + "/" + sampleName + "_Tissues_Masked_to_" + sampleObjectiveGroup + "_Transf.txt"

            print(TransformationDir)
            #sitk.ReadTransform(TransformationDir)
            sitk.ReadParameterFile(TransformationDir)
            conversionFiles.append(TransformationDir)

    nFiles = len(conversionFiles)
    print(nFiles)
    weighting = [1.0 / nFiles] * nFiles
    weightingChar = []
    for digit in weighting:
        weightingChar.append(str(digit))

    #weightingChar =  ' '.join(map(str, weighting))
    #weightingChar =  weighting
    print(weightingChar)
    print(conversionFiles)
    #p = sitk.ReadTransform(TransformationDir) #it does not work with transformation
    reference_transformation = sitk.ReadParameterFile(TransformationDir)
    p = reference_transformation#copy.deepcopy(reference_transformation)

    #strTransforms = ' '.join(conversionFiles)
    strTransforms = conversionFiles
    print(strTransforms)

    p['Transform'] = ['WeightedCombinationTransform']
    p['NormalizeCombinationWeights'] = ['true']
    p['SubTransforms'] = strTransforms

    p['TransformParameters'] = weightingChar
    p['NumberOfParameters'] = [str(nFiles)]
    p['InitialTransformParametersFileName'] = ['NoInitialTransform']
    p['ResultImagePixelType'] = ['uint8']
    p["Interpolator"] = ['NearestNeighborInterpolator']

    p["NormalizeCombinationWeights"] = ['true']

    transformixImageFilter = sitk.TransformixImageFilter()

    #sitk.PrintParameterMap(p)
    transformixImageFilter.SetTransformParameterMap(p)

    transformixImageFilter.SetMovingImage(moving_image_cropped)
    #transformixImageFilter.LogToConsoleOn()
    transformixImageFilter.LogToConsoleOff()

    resultImage = transformixImageFilter.Execute()

    print("resultImage transformation done!")

    #sitk.PrintParameterMap(p)

    transformParameterMapVector = transformixImageFilter.GetTransformParameterMap()[0]
    #sitk.PrintParameterMap(transformParameterMapVector)

    sitk.WriteParameterFile(transformParameterMapVector, TransformationDirTXT)
    sitk.WriteParameterFile(transformParameterMapVector, TransformationDirHDF5)
    sitk.WriteParameterFile(transformParameterMapVector, TransformationDirTFM)

    sitk.WriteImage(sitk.Cast(resultImage, sitk.sitkUInt8),
                    folderOutput + sampleName + '_Tissues_Masked_to_MeanTransform.tiff')

    transformAndSaveLMsElastix(moving_image_cropped, XInic, YInic, ZInic, folderOutput,
                           transformParameterMapVector, TransformationDirTXT,
                           pathLMsGMMoving, pathLMsGMModifiedForTransform, pathLMsGMCorrected, pathLMsGMMoved,
                               constDivide = constDivide)
