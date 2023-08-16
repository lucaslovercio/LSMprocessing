import cv2
import numpy as np
from PIL import Image

def functionReadTIFFMultipage(dirImage, bitdepth):
    img = Image.open(dirImage)
    images = []
    for i in range(img.n_frames):
        img.seek(i)
        images.append(np.array(img))

    del img

    height, width = np.shape(images[0])
    numImgs = len(images)
    #print(height, width, numImgs)

    if bitdepth == 8:
        volume = np.uint8(np.zeros((height, width, numImgs)))
    else:
        volume = np.uint16(np.zeros((height, width, numImgs)))

    for i in range(numImgs):
        sliceSingle = images[i]
        volume[:, :, i] = sliceSingle
    
    return volume