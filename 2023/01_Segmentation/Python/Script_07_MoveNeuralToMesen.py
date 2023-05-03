from functionReadTIFFMultipage import functionReadTIFFMultipage
from functionSaveTIFFMultipage import functionSaveTIFFMultipage
import numpy as np

#parameters
file_volume_Neural              =  "...tiff"
file_volume_Mesen               =  "...tiff"
file_volume_Tissues             =  "...tiff"

file_volume_Neural_3DSlicer     =  "...tiff"

file_volume_Neural_corrected    =  "...tiff"
file_volume_Mesen_corrected     =  "...tiff"
file_volume_Tissues_corrected   =  "...tiff"

valueNeural = 100
valueMesen = 50

#correct Neural, which is delete in Neural
volume_Neural = functionReadTIFFMultipage(file_volume_Neural, 8)
volume_Neural_3DSlicer = functionReadTIFFMultipage(file_volume_Neural_3DSlicer,8)

volume_Neural_corrected = np.zeros_like(volume_Neural)
volume_Neural_corrected = np.where(volume_Neural_3DSlicer > 0, volume_Neural_3DSlicer, 0)

functionSaveTIFFMultipage(volume_Neural_corrected, file_volume_Neural_corrected, 8)

#correct neural, which is change label from Neural to Mesen
diffVoxels = np.logical_xor(volume_Neural > 0, volume_Neural_corrected>0)

volume_Mesen_corrected = functionReadTIFFMultipage(file_volume_Mesen,8)

#extract label value
labelValue = np.max(volume_Mesen_corrected)
volume_Mesen_corrected = np.where(diffVoxels, labelValue, volume_Mesen_corrected)
functionSaveTIFFMultipage(volume_Mesen_corrected, file_volume_Mesen_corrected, 8)

#modify tissue volume. Mesen is the max value in the segmentation
volume_Tissues_corrected = functionReadTIFFMultipage(file_volume_Tissues,8)
volume_Tissues_corrected = np.where(diffVoxels, valueMesen, volume_Tissues_corrected)
functionSaveTIFFMultipage(volume_Tissues_corrected, file_volume_Tissues_corrected, 8)
