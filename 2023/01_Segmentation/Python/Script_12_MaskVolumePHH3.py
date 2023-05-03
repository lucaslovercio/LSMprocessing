import numpy as np
from functionReadTIFFMultipage import functionReadTIFFMultipage
from functionSaveTIFFMultipage import functionSaveTIFFMultipage

#parameters
file_volume_Mask =      "...tiff"
file_volume_toMask =    "...tiff"
file_dest_Masked =      "...tiff"

volumeMask = functionReadTIFFMultipage(file_volume_Mask)
volumeToMask = functionReadTIFFMultipage(file_volume_toMask)

volumeMasked = np.zeros_like(volumeToMask)
volumeMasked = np.where(volumeMask > 0, volumeToMask, volumeMasked)

functionSaveTIFFMultipage(volumeMasked, file_dest_Masked, 8)
