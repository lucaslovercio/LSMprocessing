import csv
import vtk
import SimpleITK as sitk
import numpy as np

def applyTransformToPoints(pointsVector,final_transform):
    LMs = []
    nPoints = len(pointsVector)
    for i in range(nPoints):
        pAux = pointsVector[i]
        tInv = final_transform.GetInverse()
        pAuxT = tInv.TransformPoint(pAux)
        LMs.append([pAuxT[0],pAuxT[1],pAuxT[2]])
    return LMs

def resliceVolumeITK(pathTIFFMoving, pathTIFFMovingResult, pathVolumeFixed, pathVolumeFixedCorrected,
                     landmarkTransformITK, constDivide = 1):
    spacingEmbryoCorrection = [constDivide,constDivide,constDivide]
    moving_im = sitk.ReadImage(pathTIFFMoving)
    fixed_im = sitk.ReadImage(pathVolumeFixed)

    moving_im.SetSpacing(spacingEmbryoCorrection)
    fixed_im.SetSpacing(spacingEmbryoCorrection)

    moving_im.SetOrigin(fixed_im.GetSpacing())

    itk_resampled_im = sitk.Resample(moving_im, fixed_im, landmarkTransformITK, sitk.sitkNearestNeighbor, 0)
    itk_resampled_im.SetSpacing(spacingEmbryoCorrection)
    itk_resampled_im.SetOrigin(fixed_im.GetSpacing())

    sitk.WriteImage(sitk.Cast(itk_resampled_im, sitk.sitkUInt8), pathTIFFMovingResult)


def resliceVolume(pathTIFFMoving, pathTIFFMovingResult, pathVolumeFixed, pathVolumeFixedCorrected, landmarkTransform,
                  constDivide = 0.35277777777777775 ):
    spacingEmbryoCorrection = [constDivide,constDivide,constDivide]
    #spacingEmbryoCorrection = [1, 1, 1]
    # apply to moving image
    readerVTKMoving = vtk.vtkTIFFReader()
    readerVTKMoving.SetFileName(pathTIFFMoving)
    print("Read VTK successful", pathTIFFMoving)
    #readerVTKMoving.SetOutputSpacing(spacingEmbryoCorrection[0], spacingEmbryoCorrection[1], spacingEmbryoCorrection[2])
    readerVTKMoving.Update()

    readerVTKFixed = vtk.vtkTIFFReader()
    readerVTKFixed.SetFileName(pathVolumeFixed)
    print("Read VTK successful", pathVolumeFixed)
    #readerVTKFixed.SetOutputSpacing(spacingEmbryoCorrection[0], spacingEmbryoCorrection[1], spacingEmbryoCorrection[2])
    readerVTKFixed.Update()
    #extent2 = readerVTKFixed.GetOutput().GetExtent()
    #spacingEmbryo2 = readerVTKFixed.GetOutput().GetSpacing()
    originEmbryo2 = readerVTKFixed.GetOutput().GetOrigin()

    #In the same space
    #spacingEmbryoCorrection = [0.3527, 0.3527, 0.3527]

    #spacingEmbryoCorrection = spacingEmbryo2
    VTKMoving = vtk.vtkImageReslice()
    VTKMoving.SetInputData(readerVTKMoving.GetOutput())
    # VTKMoving.SetOutputOrigin(originEmbryo2[0], originEmbryo2[1], originEmbryo2[2])
    VTKMoving.SetOutputSpacing(spacingEmbryoCorrection[0], spacingEmbryoCorrection[1], spacingEmbryoCorrection[2])
    # VTKMoving.SetOutputExtent(readerVTKFixed.GetOutput().GetExtent())
    VTKMoving.Update()

    VTKFixed = vtk.vtkImageReslice()
    VTKFixed.SetInputData(readerVTKFixed.GetOutput())
    #VTKFixed.SetOutputOrigin(originEmbryo2[0], originEmbryo2[1], originEmbryo2[2])
    VTKFixed.SetOutputSpacing(spacingEmbryoCorrection[0], spacingEmbryoCorrection[1], spacingEmbryoCorrection[2])
    VTKFixed.SetOutputExtent(readerVTKFixed.GetOutput().GetExtent())
    VTKFixed.Update()

    resliceEmbryo0 = vtk.vtkImageReslice()
    resliceEmbryo0.SetInputData(VTKMoving.GetOutput())
    resliceEmbryo0.Update()
    print("extent output: ", resliceEmbryo0.GetOutputExtent(), " SetOutputOrigin: ", resliceEmbryo0.GetOutputOrigin(),
          " GetOutputSpacing: ", resliceEmbryo0.GetOutputSpacing())

    # Applying the transformation
    mat = landmarkTransform.GetMatrix()
    transform = vtk.vtkTransform()
    transform.SetMatrix(mat)

    #resliceEmbryo0.SetOutputExtent(extent2)
    #spacingEmbryo2 = [1, 1, 1]
    resliceEmbryo0.SetOutputOrigin(originEmbryo2[0], originEmbryo2[1], originEmbryo2[2])
    resliceEmbryo0.SetOutputSpacing(spacingEmbryoCorrection[0], spacingEmbryoCorrection[1], spacingEmbryoCorrection[2])
    #resliceEmbryo0.SetOutputExtent(VTKFixed.GetOutput().GetExtent())
    resliceEmbryo0.SetOutputExtent(VTKMoving.GetOutput().GetExtent())
    resliceEmbryo0.SetInterpolationModeToNearestNeighbor()
    resliceEmbryo0.Update()


    print("extent output: ", resliceEmbryo0.GetOutputExtent(), " SetOutputSpacing: ", resliceEmbryo0.GetOutputOrigin(),
          " GetOutputSpacing: ", resliceEmbryo0.GetOutputSpacing())
    #showVolume(None, resliceEmbryo0.GetOutput(), "Before transform")
    #resliceEmbryo0.SetResliceTransform(landmarkTransform) #not working

    resliceEmbryo1 = vtk.vtkImageReslice()
    resliceEmbryo1.SetInputData(resliceEmbryo0.GetOutput())
    #resliceEmbryo1.SetInputData(VTKMoving.GetOutput())
    #resliceEmbryo1.SetResliceTransform(transform.GetInverse())
    resliceEmbryo1.SetResliceTransform(transform.GetInverse())
    resliceEmbryo1.SetOutputOrigin(originEmbryo2[0], originEmbryo2[1], originEmbryo2[2])
    resliceEmbryo1.SetOutputSpacing(spacingEmbryoCorrection[0], spacingEmbryoCorrection[1], spacingEmbryoCorrection[2])
    resliceEmbryo1.SetOutputExtent(VTKFixed.GetOutput().GetExtent())
    resliceEmbryo1.SetInterpolationModeToNearestNeighbor()
    resliceEmbryo1.Update()

    resliceEmbryo0 = []

    castFilter = vtk.vtkImageCast()
    castFilter.SetInputData(resliceEmbryo1.GetOutput())
    castFilter.SetOutputScalarTypeToUnsignedChar()
    castFilter.Update()
    print("writing Image")
    writer = vtk.vtkTIFFWriter()
    #writer.SetCompressionToLZW() #LZW compression is patented outside US so it is disabled
    writer.SetFileName(pathTIFFMovingResult)
    writer.SetInputConnection(castFilter.GetOutputPort())
    writer.Write()

    castFilter = []
    resliceEmbryo1 = []

    #Fixed image corrected
    castFilter = vtk.vtkImageCast()
    castFilter.SetInputData(VTKFixed.GetOutput())
    castFilter.SetOutputScalarTypeToUnsignedChar()
    castFilter.Update()
    print("writing Image")
    writer = vtk.vtkTIFFWriter()
    #writer.SetCompressionToLZW() #LZW compression is patented outside US so it is disabled
    writer.SetFileName(pathVolumeFixedCorrected)
    writer.SetInputConnection(castFilter.GetOutputPort())
    writer.Write()


def resample_img(itk_image, out_spacing = [0.5,0.5,0.5], is_label=True):
    #from 0.35 to 0.5
    # Resample images to 2mm spacing with SimpleITK
    original_spacing = itk_image.GetSpacing()
    original_size = itk_image.GetSize()
    print('original_spacing:' + str(original_spacing))
    print('original_size:' + str(original_size))

    out_size = [
        int(np.round(float(original_size[0]) * (original_spacing[0] / out_spacing[0]))),
        int(np.round(float(original_size[1]) * (original_spacing[1] / out_spacing[1]))),
        int(np.round(float(original_size[2]) * (original_spacing[2] / out_spacing[2])))]

    print('out_size:' + str(out_size))
    print('out_spacing:' + str(out_spacing))
    resample = sitk.ResampleImageFilter()
    resample.SetOutputSpacing(out_spacing)
    resample.SetSize(out_size)
    resample.SetOutputDirection(itk_image.GetDirection())
    resample.SetOutputOrigin(itk_image.GetOrigin())
    resample.SetTransform(sitk.Transform())
    resample.SetDefaultPixelValue(itk_image.GetPixelIDValue())

    if is_label:
        resample.SetInterpolator(sitk.sitkNearestNeighbor)
    else:
        resample.SetInterpolator(sitk.sitkBSpline)

    return resample.Execute(itk_image)