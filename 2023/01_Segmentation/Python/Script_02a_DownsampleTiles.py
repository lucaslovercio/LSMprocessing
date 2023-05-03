import os
import cv2
import numpy as np

sample = "Feb12_E115_1"
patch_size = 256
age = "E11.0"

folderOrig = "Step_01_TilesDAPI"
folderDest = "Step_02_TilesDAPI_"+str(patch_size)
ending = ".png"
filesep = '/' # / or \ depending on OS system

listing = os.listdir(folderOrig)
listing = sorted(listing)
nFiles = len(listing)
imgSize = patch_size
print("Number of tiles: " + str(nFiles))
for i in range(nFiles):
    if listing[i].endswith(ending):
        #print(listing[i])
        fullPathOrig = folderOrig + filesep + listing[i]#os.path.abspath(listing[i])
        #print(fullPathOrig)
        imgRGB = cv2.imread(fullPathOrig, cv2.IMREAD_ANYDEPTH)
        sizeRGB = np.shape(imgRGB)
        #print(sizeRGB)
        if (len(sizeRGB) >= 3):
            img = imgRGB[:,:,1]
        else:
            img = imgRGB

        img = cv2.resize(img, (imgSize, imgSize), interpolation=cv2.INTER_NEAREST)
        fullPathDest = folderDest + filesep + listing[i]
        cv2.imwrite(fullPathDest, img)


