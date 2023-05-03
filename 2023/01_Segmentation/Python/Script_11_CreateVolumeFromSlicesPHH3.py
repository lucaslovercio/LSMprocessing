import os
import numpy as np
import cv2
from PIL import Image
from PIL.TiffTags import TAGS_V2
# must also have libtiff installed in order to use lzw compression

#parameters
source_folder = ""
dest_file = "...tiff"
filesep = "/"
currentDir = os.getcwd()

#in nanometers

resX = 913.89
resY = 913.89
resZ = 4940

#processing
scale = resX/resZ

tiffs = [] #list of files (tiles) corresponding to channel
for file in os.listdir(source_folder):
    #os.chdir(currentDir + filesep + source_folder)
    os.chdir(source_folder)
    fileExtension = os.path.splitext(os.path.abspath(file))
    if fileExtension[1] == ".png":
        tiffs.append(''.join(fileExtension))
tiffs = sorted(tiffs)

numImgs = len(tiffs) #number of tiles

multifile_tiff = []
for i in range(0, numImgs):
    full_path_source = tiffs[i]
    slice_original = cv2.imread(full_path_source, cv2.IMREAD_GRAYSCALE)
    
    height, width = np.shape(slice_original)
    resized = (int(scale * width), int(scale * height))
    img_downsample = cv2.resize(slice_original, resized, cv2.INTER_NEAREST)
    img_downsample = np.uint8(img_downsample)
    img_downsample = Image.fromarray(img_downsample)
    multifile_tiff.append(img_downsample)


os.chdir(currentDir)
firstSlice = multifile_tiff[0]
firstSlice.save(dest_file, save_all = True, append_images = multifile_tiff[1:], compression = "tiff_lzw")

print("Finished!")

