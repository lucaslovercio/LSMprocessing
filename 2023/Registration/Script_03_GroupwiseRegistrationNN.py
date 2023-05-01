import os
import SimpleITK as sitk
#import gui
import registration_gui as rgui
from functionsRegistration import resample_img

#import vtk
#from functionsRegistration import convertLMsToVTKPoints, resliceVolume, saveLMs, convertVTKPointsToLMs, applyTransformToPoints, modifyOrigin

#------------------INPUTS-------------------

workFolder = "Registration/"
group = "E11.0_1"
Downsampled = False
sampleObjectiveGroup= "Feb12_E115_2"

RegType = "AffineGroupwiseNN"
previousReg = "LM38RigidReg_itk"

#PIL TIFFs with dimension 1, but MATLAB generate TIFFs with 0.35...
#Set according to the file. Landmarks are always using voxel size = 1
constDivide = 1

sampleNames = ["Feb3_E11_4_v2","Feb3_E11_4_v2_flipped","Feb12_E115_1","Feb12_E115_1_flipped",\
               "Feb12_E115_2","Feb12_E115_2_flipped","Feb12_E115_8","Feb12_E115_8_flipped",\
               "May10_E11_2","May10_E11_2_flipped"]

for sampleObjectiveGroup in sampleNames:
    print("----------Moving-----------------------------")
    print(sampleObjectiveGroup)

    # --------------OBJECTIVE-------------------------

    pathLMsNeckObjective = workFolder + group + "/" + previousReg + "/" + sampleObjectiveGroup + "_to_Objective" + "/" \
                           + sampleObjectiveGroup + "_Tissues_Masked_neck_Objective.pts"
    pathLMsFaceObjective = workFolder + group + "/" + previousReg + "/" + sampleObjectiveGroup + "_to_Objective" + "/" \
                           + sampleObjectiveGroup + "_Tissues_Masked_face_Objective.pts"
    pathLMsGMObjective = workFolder + group + "/" + previousReg + "/" + sampleObjectiveGroup + "_to_Objective" + "/" \
                         + sampleObjectiveGroup + "_Tissues_Masked_38lms_Objective.pts"
    pathVolumeObjective = workFolder + group + "/" + previousReg + "/" + sampleObjectiveGroup + "_to_Objective" + "/" \
                          + sampleObjectiveGroup + "_Tissues_Reg.tiff"
    pathVolumeObjectiveCropped = workFolder + group + "/" + previousReg + "/" + sampleObjectiveGroup + "_to_Objective" + "/" \
                                 + sampleObjectiveGroup + "_Tissues_Objective_Cropped.tiff"

    if Downsampled:

        pathVolumeObjective = workFolder + group + "/" + previousReg + "/" + "Downsampled" + "/" + sampleObjectiveGroup + "_to_Objective" + "/" \
                              + sampleObjectiveGroup + "_downsampled_Tissues_Reg.tiff"
        pathVolumeObjectiveCropped = workFolder + group + "/" + previousReg + "/" + "Downsampled" + "/" + sampleObjectiveGroup + "_to_Objective" + "/" \
                                     + sampleObjectiveGroup + "_downsampled_Tissues_Objective_Cropped.tiff"


    print(pathVolumeObjective)

    for sampleName in sampleNames:
        if (1 == 1):
            # if not(sampleName == sampleObjectiveGroup):
            #sampleName =    "Dec2_E10_16_flipped"
            print("----------Moving-----------------------------")
            print(sampleName)

            #-------------- Moving ------------------------------

            folder =        workFolder+ group +"/"+ previousReg + "/"+sampleName+"_to_Objective/"
            regTypeFolder = workFolder+ group +"/"+ previousReg + "/"+ RegType+"/"

            if Downsampled:
                folder = workFolder + group + "/" + previousReg + "/" + "Downsampled" + "/" + sampleName + "_to_Objective/"
                regTypeFolder = workFolder + group + "/" + previousReg  + "/" + "Downsampled" + "/" + RegType + "/"

            try:
                os.mkdir(regTypeFolder)
            except OSError as error:
                print(error)

            folderOutput =  regTypeFolder +sampleName+"_to_"+ sampleObjectiveGroup

            try:
                os.mkdir(folderOutput)
            except OSError as error:
                print(error)

            pathLMsFaceMoving = folder + sampleName+"_Tissues_Masked_face_to_"+sampleObjectiveGroup+".pts"
            pathLMsNeckMoving = folder + sampleName+"_Tissues_Masked_neck_to_"+sampleObjectiveGroup+".pts"
            pathLMsGMMoving =   folder + sampleName+"_Tissues_Masked_38lms_to_"+sampleObjectiveGroup+".pts"

            pathLMsFaceCorrected = folderOutput  + "/"+ sampleName + "_Tissues_Masked_face_Corrected.pts"
            pathLMsNeckCorrected = folderOutput  + "/"+ sampleName + "_Tissues_Masked_neck_Corrected.pts"
            pathLMsGMCorrected = folderOutput  + "/"+ sampleName + "_Tissues_Masked_38lms_Corrected.pts"

            pathLMsFaceMoved = folderOutput  + "/"+ sampleName+"_Tissues_Masked_face_to_"+sampleObjectiveGroup+".pts"
            pathLMsNeckMoved = folderOutput  + "/"+ sampleName+"_Tissues_Masked_neck_to_"+sampleObjectiveGroup+".pts"
            pathLMsGMMoved = folderOutput  + "/"+ sampleName+"_Tissues_Masked_38lms_to_"+sampleObjectiveGroup+".pts"

            #-------------- Moving ------------------------------

            pathVolumeFixedCorrected = folderOutput + "/"+"FixedCorrected.tiff"

            MovingImage = folder + sampleName + "_Tissues_Reg.tiff"

            if Downsampled:
                MovingImage = folder + sampleName + "_downsampled_Tissues_Reg.tiff"

            TransformationDirTXT = folderOutput  + "/"+ sampleName + "_Tissues_Masked_to_"+sampleObjectiveGroup+"_Transf.txt"
            TransformationDirHDF5 = folderOutput + "/" + sampleName + "_Tissues_Masked_to_" + sampleObjectiveGroup + "_Transf.hdf5"
            TransformationDirTFM = folderOutput + "/" + sampleName + "_Tissues_Masked_to_" + sampleObjectiveGroup + "_Transf.tfm"
            LogFileName =           folderOutput + "/" + sampleName + "_Tissues_Masked_to_" + sampleObjectiveGroup + "_Log.txt"

            #------------------------------------------------

            # Intensity Registration

            fixed_image = sitk.ReadImage(pathVolumeObjective, sitk.sitkFloat32)
            print("-------Fixed image stats------------")
            #Crop the image because it is too big
            print(fixed_image.GetOrigin())
            imageSize = fixed_image.GetSize()

            #Between E10.0 and E11.0
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


            fixed_image.SetSpacing([constDivide,constDivide,constDivide])

            fixed_image_cropped = fixed_image[XInic:(imageSize[0]-XInic), YInic:imageSize[1], ZInic:(imageSize[2]-ZInic)]

            #print("-------Fixed Image Cropped------------")
            del fixed_image

            #if group == "E11.5":
                #Resize to fit in memory
                #fixed_image_cropped = resample_img(fixed_image_cropped)


            sitk.WriteImage(sitk.Cast(fixed_image_cropped, sitk.sitkUInt8),pathVolumeObjectiveCropped)
            #print("saving: " + pathVolumeObjectiveCropped)

            moving_image = sitk.ReadImage(MovingImage, sitk.sitkFloat32)

            moving_image.SetSpacing([constDivide, constDivide, constDivide])
            #print(moving_image.GetSpacing())

            moving_image_cropped = moving_image[XInic:(imageSize[0]-XInic), YInic:imageSize[1], ZInic:(imageSize[2]-ZInic)]
            #print("-------After Image Cropped------------")

            del moving_image

            #if group == "E11.5":
                #Resize to fit in memory
                #moving_image_cropped = resample_img(moving_image_cropped)

            pathVolumeMovingCropped = folderOutput + "/"+ sampleName + "_Tissues_Masked_Cropped.tiff"
            #sitk.WriteImage(sitk.Cast(moving_image_cropped, sitk.sitkUInt8), pathVolumeMovingCropped) #SimpleElastix did not implemented compression
            #print("saving: " + pathVolumeMovingCropped)

            #-----Registration-------
            #sitk.SetParameter("WriteIterationInfo", ["true"])

            elastixImageFilter = sitk.ElastixImageFilter()
            elastixImageFilter.SetParameter("WriteIterationInfo", ["true"])
            elastixImageFilter.LogToConsoleOff()
            elastixImageFilter.LogToFileOn()
            #elastixImageFilter.SetLogFileName(LogFileName)
            #elastixImageFilter.SetLogToFile(true)
            elastixImageFilter.SetOutputDirectory(folderOutput) #Elastix Logger has bugs https://github.com/SuperElastix/SimpleElastix/issues/267. It is not working

            elastixImageFilter.SetFixedImage(fixed_image_cropped)
            elastixImageFilter.SetMovingImage(moving_image_cropped)

            parameterMap = sitk.GetDefaultParameterMap("affine")
            #sitk.PrintParameterMap(parameterMap)
            #parameterMap["Metric"] = ['AdvancedMeanSquares'] #be careful with the
            parameterMap['MaximumNumberOfIterations'] = ['30000'] #30000
            parameterMap["Interpolator"] = ['LinearInterpolator'] #Interpolators are used to interpolate intensities of images and deformations of transforms.
            parameterMap["ResampleInterpolator"] = ["FinalNearestNeighborInterpolator"] #ResampleInterpolators are used to interpolate the intensity of an image when resampling it.
            parameterMap["NumberOfResolutions"] = ['5']
            parameterMap["FixedImagePyramid"] = ['FixedSmoothingImagePyramid']
            parameterMap["Registration"] = ['MultiResolutionRegistration']
            #sitk.PrintParameterMap(parameterMap)

            elastixImageFilter.SetParameterMap(parameterMap)

            elastixImageFilter.Execute()

            transformParameterMapVector = elastixImageFilter.GetTransformParameterMap()[0]
            #sitk.PrintParameterMap(transformParameterMapVector)

            sitk.WriteParameterFile(transformParameterMapVector, TransformationDirTXT)
            sitk.WriteParameterFile(transformParameterMapVector, TransformationDirHDF5)
            sitk.WriteParameterFile(transformParameterMapVector, TransformationDirTFM)
            #sitk.WriteTransform()

            del fixed_image_cropped

            #To be sure that it does not modify intensities
            transformParameterMapVector["Interpolator"] = ['NearestNeighborInterpolator']
            transformParameterMapVector["ResampleInterpolator"] = ["FinalNearestNeighborInterpolator"]

            print("-------------Entering test-------------")
            resultImage = sitk.Transformix(moving_image_cropped, transformParameterMapVector)

            del moving_image_cropped

            #Write resulting moved image to verify. Not recommended to use regularly, just for debugging.
            #sitk.WriteImage(sitk.Cast(elastixImageFilter.GetResultImage(), sitk.sitkUInt8),folderOutput  + "/"+ sampleName + "_Tissues_Masked_to_"+sampleObjectiveGroup+'.tiff')

            del elastixImageFilter
            del folderOutput





