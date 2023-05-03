from skimage import exposure
from skimage import io
from skimage import img_as_ubyte
import os
import cv2
import numpy as np
from functionPercNorm import functionPercNorm

folderVolume =      '' # folder containing slices
folderDestTiles =   'Step_01_TilesDAPI' # destination folder
sampleName =        'Feb12_E115_7' # prefix for the resulting files
filesep = '/' # / or \ depending on OS system
outputName = os.path.abspath(folderDestTiles) + filesep + sampleName 

# output parameters
# nbits = 8, can be 16
nChannels = 1 # can be 3
fileFormat = '.png'

patchSize = 1024
spx = np.uint16(patchSize * 0.1)

#------------------------------------------------------------

listing = os.listdir(folderVolume)
listing = sorted(listing)
nFiles = len(listing)
usefulSize = patchSize - (2*spx)
currentDir = os.getcwd()

for imgNumber in range(nFiles):
    os.chdir(folderVolume)
    sliceNumber = imgNumber + 1
    print('Processing:' + str(imgNumber + 1))
    fullpathOrig = os.path.abspath(listing[imgNumber]) 
    print(cv2.haveImageReader(fullpathOrig))
    img = cv2.imread(fullpathOrig, cv2.IMREAD_ANYDEPTH)
    hOrig, wOrig = np.shape(img)
    #
    imgNorm = functionPercNorm(img) #Normalization fixed
    imgNorm = exposure.rescale_intensity(imgNorm, out_range=(0, 255))
    #
    hImgAugmented = hOrig + patchSize
    wImgAugmented = wOrig + patchSize
    imgAugmented = np.zeros([hImgAugmented, wImgAugmented])
    #
    imgAugmented[(spx + 1):(spx + hOrig + 1), (spx +1):(spx + wOrig + 1)] = imgNorm
    
    # Saving the tiles
    numBlocksR = 0
    for r in range(0, hOrig , usefulSize): #iteration inside the true image
        numBlocksC = 0
        numBlocksR += 1
        for c in range(0, wOrig, usefulSize):
            numBlocksC += 1
            #add neighbourhood to tile

            # top left corner
            r1 = r
            c1 = c

            #upper neighbourhood in the true image
            r0 = r1 - spx
            c0 = c1 - spx

            #bottom right neighbourhood in the true image
            r2 = r0 + patchSize - 1
            c2 = c0 + patchSize - 1

            #translation to augmented image
            r0_a = r0 + spx
            c0_a = c0 + spx
            r2_a = r2 + spx
            c2_a = c2 + spx

            block = imgAugmented[r0_a:(r2_a + 1), c0_a:(c2_a + 1)]

            name = (outputName + '_slice_' + str(sliceNumber).zfill(4) + '_block_' +
                str(numBlocksR).zfill(2) + '_' + str(numBlocksC).zfill(2) + fileFormat)
            if nChannels == 3:
                rgb = np.concatenate(3, block, block, block, block)
            else:
                rgb = block

            os.chdir(folderDestTiles)
            cv2.imwrite(name, rgb.astype(np.uint8)) #saves as 8 bit image
            #io.imsave(name, img_as_ubyte(rgb), check_contrast = False)
              
    # save slice data
    infoName = outputName + '_slice_' + str(sliceNumber).zfill(4) + '_info.txt'

    fid = open(infoName , 'w')
    fid.write('sampleName:'+ sampleName + 
        '\nimgNumber:' + str(sliceNumber) +
        '\nhOrig:' + str(hOrig) +
        ' \nwOrig:' + str(wOrig) +
        '\nhImgAugmented:' + str(hImgAugmented) +
        '\nwImgAugmented:' + str(wImgAugmented) +
        '\npatchSize:' + str(patchSize) +
        '\nnumBlocksR:' + str(numBlocksR) +
        '\nnumBlocksC:' + str(numBlocksC) +
        '\nspx:' + str(spx) +
        '\nusefulSize:' + str(usefulSize))
    
    fid.close()
    
