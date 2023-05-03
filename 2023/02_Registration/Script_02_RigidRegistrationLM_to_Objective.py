import os
import SimpleITK as sitk
from functionsRegistration import resliceVolumeITK
import functionsLandmarks as flms

#------------------INPUTS-------------------

workFolder = "Registration/"
group = "E11.0_1"
sampleObjectiveGroup= "Feb12_E115_2"

#PIL TIFFs with dimension 1, but MATLAB generate TIFFs with 0.35...
#Set according to the file. Landmarks are always using voxel size = 1
constDivide = 1

RegType = "LM38RigidReg_itk"

sampleNames = ["Feb3_E11_4_v2","Feb3_E11_4_v2_flipped","Feb12_E115_1","Feb12_E115_1_flipped",\
               "Feb12_E115_2","Feb12_E115_2_flipped","Feb12_E115_8","Feb12_E115_8_flipped",\
               "May10_E11_2","May10_E11_2_flipped"]

for sampleName in sampleNames:
    #sampleName =    "Dec2_E10_16_flipped"
    print("---------------------------------------")
    print(sampleName)

    #--------------OBJECTIVE-------------------------

    pathLMsNeckObjective = workFolder + group + "/" + sampleObjectiveGroup + "_Objective" + "/" \
        + sampleObjectiveGroup + "_Tissues_Masked_neck_Objective.csv"
    pathLMsFaceObjective = workFolder + group + "/" + sampleObjectiveGroup + "_Objective" + "/" \
        + sampleObjectiveGroup + "_Tissues_Masked_face_Objective.csv"
    pathLMsGMObjective = workFolder + group + "/" + sampleObjectiveGroup + "_Objective" + "/" \
        + sampleObjectiveGroup + "_Tissues_Masked_38lms_Objective.csv"
    pathVolumeObjective = workFolder + group + "/" + sampleObjectiveGroup + "_Objective" + "/" \
        + sampleObjectiveGroup + "_Tissues_Objective.tiff"

    #-------------- Moving ------------------------------

    folder =        workFolder+ group +"/"+sampleName+"/"
    regTypeFolder = workFolder+ group +"/"+ RegType+"/"

    try:
        os.mkdir(regTypeFolder)
    except OSError as error:
        print(error)

    folderOutput =  regTypeFolder +sampleName+"_to_Objective/"

    try:
        os.mkdir(folderOutput)
    except OSError as error:
        print(error)

    pathLMsFaceMoving = folder + sampleName+"_Tissues_Masked_face.csv"
    pathLMsNeckMoving = folder + sampleName+"_Tissues_Masked_neck.csv"
    pathLMsGMMoving =   folder + sampleName+"_Tissues_Masked_38lms.csv"

    pathLMsFaceMoved = folderOutput + sampleName+"_Tissues_Masked_face_to_Objective.csv"
    pathLMsNeckMoved = folderOutput + sampleName+"_Tissues_Masked_neck_to_Objective.csv"
    pathLMsGMMoved = folderOutput + sampleName+"_Tissues_Masked_38lms_to_Objective.csv"

    pathTransformationTXT = folderOutput + sampleName + "_transf_to_Objective.txt"
    pathTransformationTFM = folderOutput + sampleName + "_transf_to_Objective.tfm"

    pathTransformationInvTXT = folderOutput + sampleName + "_transf_inverse_to_Objective.txt"
    pathTransformationInvTFM = folderOutput + sampleName + "_transf_inverse_to_Objective.tfm"

    #-------------- Moving ------------------------------

    pathVolumeFixedCorrected = folderOutput+"FixedCorrected.tiff"

    pathMovingVolumes = [\
                        folder + sampleName + "_Tissues.tiff", \
                        folder + sampleName + "_ProliferatingCells_Masked.tiff", \
                        #folder + sampleName + "_Neural.tiff", \
                        #folder + sampleName + "_Mesen.tiff", \
                        ]

    pathMovedVolumes = [\
                        folderOutput + sampleName + "_Tissues_Reg.tiff", \
                        folderOutput + sampleName + "_ProliferatingCells_Masked_Reg.tiff", \
                        #folderOutput + sampleName + "_Neural_Reg.tiff", \
                        #folderOutput + sampleName + "_Mesen_Reg.tiff", \
                        ]

    #------------------------------------------------

    #Get spacing from volumes
    #Read volumes

    moving_im = sitk.ReadImage(pathMovingVolumes[0])
    spacingEmbryoMoving = moving_im.GetSpacing()
    print("spacingEmbryoMoving" + str(spacingEmbryoMoving))

    fixed_im = sitk.ReadImage(pathVolumeObjective)
    print("Read ITK successful", pathVolumeObjective)

    extent2 = fixed_im
    spacingEmbryoFixed = fixed_im.GetSpacing()
    print("spacingEmbryoFixed" + str(spacingEmbryoFixed))
    #Registration of landmarks to the fixed

    sourcePoints = flms.convertLMsToITKPoints(pathLMsGMMoving, constDivide = constDivide)
    targetPoints = flms.convertLMsToITKPoints(pathLMsGMObjective , constDivide = constDivide)

    #For E9.5, we do not have GM LMs, we use the Face LMS (five)
    # sourcePoints = flms.convertLMsToITKPoints(pathLMsFaceMoving)
    # targetPoints = flms.convertLMsToITKPoints(pathLMsFaceObjective)

    #landmarkTransformITK = sitk.LandmarkBasedTransformInitializer(sitk.VersorRigid3DTransform(),
    #                                                              targetPoints,
    #                                                              sourcePoints,fixed_im)
    InitTx = sitk.VersorRigid3DTransform()
    landmarkTransformITKFilter = sitk.LandmarkBasedTransformInitializerFilter()
    landmarkTransformITKFilter.SetFixedLandmarks(targetPoints)
    landmarkTransformITKFilter.SetMovingLandmarks(sourcePoints)
    landmarkTransformITKFilter.SetReferenceImage(fixed_im)
    landmarkTransformITK = landmarkTransformITKFilter.Execute(InitTx)

    sitk.WriteTransform(landmarkTransformITK, pathTransformationTXT)
    sitk.WriteTransform(landmarkTransformITK, pathTransformationTFM)

    # Display the transformation matrix that was computed
    #mat = landmarkTransformITK.GetMatrix()
    #print(mat)

    #For landmarks, it is neccesary to apply the inverse of the transformation! Even if it is computed with landmarks. Images match!
    landmarkTransformITKInverse = landmarkTransformITK.GetInverse()
    sitk.WriteTransform(landmarkTransformITKInverse, pathTransformationInvTXT)
    sitk.WriteTransform(landmarkTransformITKInverse, pathTransformationInvTFM)

    #convert and save LMs moved
    #Neck landmarks
    if os.path.exists(pathLMsNeckMoving):
        print('-----Saving Neck LMs Transformed------')
        neckPoints = flms.getLMs(pathLMsNeckMoving)
        LMsNeckMoved = [landmarkTransformITKInverse.TransformPoint(p) for p in neckPoints]
        flms.saveLMs(pathLMsNeckMoved, LMsNeckMoved)
    else:
        print('-----No Neck LMs------')

    #Face landmarks
    if os.path.exists(pathLMsFaceMoving):
        print('-----Saving Face LMs Transformed------')
        facePoints = flms.getLMs(pathLMsFaceMoving)
        LMsFaceMoved = [landmarkTransformITKInverse.TransformPoint(p) for p in facePoints]
        flms.saveLMs(pathLMsFaceMoved, LMsFaceMoved)
    else:
        print('-----No Face LMs------')

    #GM landmarks
    if os.path.exists(pathLMsGMMoving):
        print('-----Saving 38 LMs Transformed------')
        GMPoints = flms.getLMs(pathLMsGMMoving, constDivide = constDivide)
        LMsGMMoved = [landmarkTransformITKInverse.TransformPoint(p) for p in GMPoints]
        flms.saveLMs(pathLMsGMMoved, LMsGMMoved, constDivide = constDivide)
    else:
        print('-----No 38 LMs------')

    nVolumes = len(pathMovingVolumes)

    for i in range(nVolumes):
        pathMovingVolume = pathMovingVolumes[i]
        pathResultVolume = pathMovedVolumes[i]
        resliceVolumeITK(pathMovingVolume, pathResultVolume, pathVolumeObjective, pathVolumeFixedCorrected,
                         landmarkTransformITK, constDivide = constDivide)


