import os
import numpy as np
import cv2
from PIL import Image
#parameters

sample = "Feb12_E115_1"
patch_size = "1024"

source_folder_segmented =""
dest_file_neural = 	 "...tiff"
dest_file_mesen = 	 "...tiff"

currentDir = os.getcwd()
filesep = "/"

valueNeural = 100
valueMesen = 50

#in nanometers
#search the voxel dimensions in Superhero Embryos, Sheet 2
resX = 913.89
resY = 913.89
resZ = 4940

#----------------------------------------------------------------

#processing
scale = resX/resZ

tiffs = [] #list of files (tiles) corresponding to channel
for file in os.listdir(source_folder_segmented):
    #os.chdir(currentDir + filesep + source_folder_segmented)
    os.chdir(source_folder_segmented)
    fileExtension = os.path.splitext(os.path.abspath(file))
    if fileExtension[1] == ".png":
        tiffs.append(''.join(fileExtension))
tiffs = sorted(tiffs)

numImgs = len(tiffs) #number of tiles

messenTiffs = []
neuralTiffs = []
for i in range(numImgs):
    full_path_source = tiffs[i]
    slice_original = cv2.imread(full_path_source, cv2.IMREAD_GRAYSCALE)

    #tissue thresholding
    #sliceNeural = slice_original == valueNeural
    #sliceMesen = slice_original == valueMesen

    sliceNeural = np.zeros_like(slice_original)
    sliceNeural = np.where(slice_original == valueNeural, 1, 0)
    sliceMesen = np.zeros_like(slice_original)
    sliceMesen = np.where(slice_original == valueMesen, 1, 0)
    
    height, width = np.shape(slice_original)
    resized = [int(scale * width), int(scale * height)]
  
    sliceMesen = Image.fromarray(np.uint8(sliceMesen))
    img_downsample_mesen = sliceMesen.resize(resized, resample=Image.NEAREST)
    messenTiffs.append(img_downsample_mesen)

    sliceNeural = Image.fromarray(np.uint8(sliceNeural))
    img_downsample_neural = sliceNeural.resize(resized, resample=Image.NEAREST)
    neuralTiffs.append(img_downsample_neural)

os.chdir(currentDir)
messenTiffs[0].save(dest_file_mesen, save_all = True, append_images = messenTiffs[1:], compression = "tiff_lzw")
neuralTiffs[0].save(dest_file_neural, save_all = True, append_images = neuralTiffs[1:], compression = "tiff_lzw")

print("Finished!")
