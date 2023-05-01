import os
import vtk
from functionsRegistration import resliceVolume
from functionsLandmarks import convertLMsToVTKPoints, saveLMs, convertVTKPointsToLMs

group = "E11.0_1"

sampleName= "Feb12_E115_2"
folder =        "Registration/"+group+"/"+sampleName+"/"
folderOutput =  "Registration/"+group+"/"+sampleName+"_Objective/"

#PIL TIFFs with dimension 1, but MATLAB generate TIFFs with 0.35...
#Set according to the file. Landmarks are always using voxel size = 1
constDivide = 1

try:
    os.mkdir(folderOutput)
except OSError as error:
    print(error)

pathLMsFaceMoving =                    folder + sampleName+"_Tissues_Masked_face.csv"
pathLMsNeckMoving =                    folder + sampleName+"_Tissues_Masked_neck.csv"
pathLMsGMMoving =                    folder + sampleName+"_Tissues_Masked_38lms.csv"

pathLMsNeck_to_Reference = folderOutput + sampleName+"_Tissues_Masked_neck_Objective.csv"
pathLMsFace_to_Reference = folderOutput + sampleName+"_Tissues_Masked_face_Objective.csv"
pathLMsGM_to_Reference = folderOutput + sampleName+"_Tissues_Masked_38lms_Objective.csv"

pathLMsNeckFixed =     "ReferenceVolume/Reference_volume_neck_base.csv"
pathVolumeFixed =      "ReferenceVolume/Reference_volume_neck_base.tiff"
pathVolumeFixedCorrected = folderOutput+"FixedCorrected.tiff"

pathMovingVolumes = [folder + sampleName + "_Tissues.tiff", \
                     folder + sampleName + "_Neural.tiff", \
                     folder + sampleName + "_Mesen.tiff", \
                     ]

pathMovedVolumes = [
                    folderOutput + sampleName + "_Tissues_Objective.tiff", \
                    folderOutput + sampleName + "_Neural.tiff", \
                    folderOutput + sampleName + "_Mesen.tiff", \
                    ]

#------------------------------------------------

#Get spacing from volumes
#Read volumes
readerVTKMoving = vtk.vtkTIFFReader()
readerVTKMoving.SetFileName(pathMovingVolumes[0])
print("Read VTK successful", pathMovingVolumes[0])
readerVTKMoving.Update()
spacingEmbryoMoving = readerVTKMoving.GetOutput().GetSpacing()
readerVTKMoving = []
print("spacingEmbryoMoving" + str(spacingEmbryoMoving))

readerVTKFixed = vtk.vtkTIFFReader()
readerVTKFixed.SetFileName(pathVolumeFixed)
print("Read VTK successful", pathVolumeFixed)
readerVTKFixed.Update()
extent2 = readerVTKFixed.GetOutput().GetExtent()
spacingEmbryoFixed = readerVTKFixed.GetOutput().GetSpacing()
readerVTKFixed = []
print("spacingEmbryoFixed" + str(spacingEmbryoFixed))
#Registration of landmarks to the fixed

sourcePoints = convertLMsToVTKPoints(pathLMsNeckMoving, constDivide = constDivide)
targetPoints = convertLMsToVTKPoints(pathLMsNeckFixed, constDivide = constDivide)


landmarkTransformVTK = vtk.vtkLandmarkTransform()
landmarkTransformVTK.SetSourceLandmarks(sourcePoints)
landmarkTransformVTK.SetTargetLandmarks(targetPoints)
landmarkTransformVTK.SetModeToRigidBody()
landmarkTransformVTK.Update()

# Display the transformation matrix that was computed
mat = landmarkTransformVTK.GetMatrix()
print(mat)

#convert and save LMs moved
#Neck landmarks
print('-----Saving Neck LMs Transformed------')
lmsNeck_to_Reference = vtk.vtkPoints()
landmarkTransformVTK.TransformPoints(sourcePoints,lmsNeck_to_Reference)
LMsNeckMoved = convertVTKPointsToLMs(lmsNeck_to_Reference, constDivide = constDivide)
saveLMs(pathLMsNeck_to_Reference,LMsNeckMoved, constDivide = constDivide)

#Face landmarks
print('-----Saving Face LMs Transformed------')
lmsFace_to_Reference = vtk.vtkPoints()
facePoints = convertLMsToVTKPoints(pathLMsFaceMoving)
landmarkTransformVTK.TransformPoints(facePoints,lmsFace_to_Reference)
LMsFaceMoved = convertVTKPointsToLMs(lmsFace_to_Reference, constDivide = constDivide)
saveLMs(pathLMsFace_to_Reference,LMsFaceMoved, constDivide = constDivide)

lmsGM_to_Reference = vtk.vtkPoints()
GMPoints = convertLMsToVTKPoints(pathLMsGMMoving)
landmarkTransformVTK.TransformPoints(GMPoints,lmsGM_to_Reference)
LMsGMMoved = convertVTKPointsToLMs(lmsGM_to_Reference, constDivide = constDivide)
saveLMs(pathLMsGM_to_Reference,LMsGMMoved, constDivide = constDivide)

nVolumes = len(pathMovingVolumes)

for i in range(nVolumes):
    pathMovingVolume = pathMovingVolumes[i]
    pathResultVolume = pathMovedVolumes[i]
    resliceVolume(pathMovingVolume, pathResultVolume, pathVolumeFixed, pathVolumeFixedCorrected, landmarkTransformVTK,
                  constDivide = constDivide)


