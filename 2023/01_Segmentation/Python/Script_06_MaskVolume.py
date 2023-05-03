import numpy as np
from functionReadTIFFMultipage import functionReadTIFFMultipage
from functionSaveTIFFMultipage import functionSaveTIFFMultipage

#parameters
file_volume_Mask =  "...tiff"
file_volume_toMask ="...tiff"
file_dest_Masked =  "...tiff"

volumeMask = functionReadTIFFMultipage(file_volume_Mask, 8)
volumeToMask = functionReadTIFFMultipage(file_volume_toMask, 8)

volumeMasked = np.zeros_like(volumeToMask)
volumeMasked = np.where(volumeMask > 0, volumeToMask, volumeMasked)

functionSaveTIFFMultipage(volumeMasked, file_dest_Masked, 8)
