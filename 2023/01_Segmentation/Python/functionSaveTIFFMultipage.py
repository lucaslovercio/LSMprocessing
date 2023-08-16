import numpy as np
from PIL import Image

def functionSaveTIFFMultipage(volume, fileNameOutput, bitdepth):
    h, w, d = np.shape(volume)

    if bitdepth == 8:
        volume = np.uint8(volume)
    else:
        volume = np.uint16(volume)

    imlist = []
    #for m in volume:
    for i in range(d):
        #imlist.append(Image.fromarray(m))
        m = volume[:,:,i]
        imlist.append(Image.fromarray(m))

    imlist[0].save(fileNameOutput, compression="tiff_lzw", save_all=True,
                   append_images=imlist[1:])
