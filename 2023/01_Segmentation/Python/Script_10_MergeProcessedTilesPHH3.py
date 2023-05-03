import numpy as np
import os
import cv2

folderTilesOriginal =   "Step_08_TilesPHH3" # folder with unprocessed tiles
folderTilesProcessed =  "Step_09_TilesPHH3_Processed" # folder with segmented tiles using any U-net
folderOutputMerged =    "Step_10_ProliferationSlices"  # folder where the restored segmented slices will be saved
filesep = '/' # / or \ depending on OS system
mainDir = os.getcwd()

#output parameters
nbits = 8 # can be 16
fileformat = ".png"
#------------------------------------------------

# look for all the .txt where information is stored
filesToSearch = []
for file in os.listdir(folderTilesOriginal):
    fileExtension = os.path.splitext(os.path.abspath(file))
    if fileExtension[1] == ".txt":
        filesToSearch.append(''.join(fileExtension))
nFiles = len(filesToSearch)
filesToSearch = sorted(filesToSearch)

# for each file in the directory
for i in range(nFiles):
    os.chdir(folderTilesOriginal)
    filename = os.path.basename(filesToSearch[i])
    #print(filename)
    with open(filename) as file:
        lines = [line.rstrip() for line in file] 
    file.close()
    lines = [line.split(':', 1) for line in lines]
    lines = np.array(lines).reshape(11, 2)
    
    #read the ouput name and imgNumber
    sampleName = lines[0][1]
    #print(sampleName)
    sliceNumber = int(lines[1][1])
    #print(sliceNumber)
    hOrig = int(lines[2][1])
    #print(hOrig)
    #print("hOrig:" + str(hOrig))
    wOrig = int(lines[3][1])
    #print("wOrig:" + str(wOrig))
    hImgAugmented = int(lines[4][1])
    wImgAugmented = int(lines[5][1])
    patchSize = int(lines[6][1])
    #print("patchSize:"+str(patchSize))
    numBlocksR = int(lines[7][1])
    #print("numBlocksR:"+str(numBlocksR))
    numBlocksC = int(lines[8][1])
    #print("numBlocksC:" + str(numBlocksC))
    spx = int(lines[9][1])
    #print("spx:" + str(spx))
    usefulSize = int(lines[10][1])
    #print("usefulSize:" + str(usefulSize))

    os.chdir(folderTilesProcessed)
    # build the first part of the string
    firstPatPath = os.getcwd() + filesep + sampleName + "_slice_" + str(sliceNumber).zfill(4)
    #print(firstPatPath)
    print("Processing: " + str(sliceNumber).zfill(4))

    # read the original hOrig and wOrig, create an empty matirx of zeros
    if nbits == 8:
        mergedImg = np.uint8(np.zeros([hOrig, wOrig]))
    else:
        mergedImg = np.uint16(np.zeros([hOrig, wOrig]))
    
    # read the numBlocksR and numBlocksC
    for r in range(1, numBlocksR+1):
        #print("row:"+str(r))
       for c in range(1, numBlocksC+1):
            #print("col:" + str(c))
            # for r in Rows and c in Cols
            # reconstruct te png filename
            processedFile = firstPatPath + "_block_" + str(r).zfill(2) + "_" + str(c).zfill(2) + "_cp_masks" + fileformat
            #print(processedFile)
            # read the file
            img = cv2.imread(processedFile, cv2.IMREAD_ANYDEPTH)
            #print(np.shape(img))
            hProcessed, wProcessed = np.shape(img)
            if hProcessed != patchSize:
                img = cv2.resize(img, [patchSize, patchSize], interpolation = cv2.INTER_NEAREST)

            #img = img[(spx + 1):(len(img) - spx + 1), (spx + 1):(len(img[0]) - spx + 1)]

            img = img[spx:(patchSize - spx), spx:(patchSize - spx)]
            #print("after cropping overlapping:" + str(np.shape(img)))
            if r ==1 and c==1: # first iteration
                if img.dtype == np.uint16:
                    mergedImg = mergedImg.astype(np.uint16)
                else:
                    mergedImg = mergedImg.astype(np.uint8) # only one 8bits channel
            
            # multiply the r and c per patchSize
            rOriginal = (r-1) * usefulSize
            cOriginal = (c-1) * usefulSize

            #print("inic rOriginal:" + str(rOriginal) + " cOriginal:" + str(cOriginal) )

            endRowOriginal = rOriginal + usefulSize
            endColOriginal = cOriginal + usefulSize
            endRowPatch = usefulSize
            endColPatch = usefulSize

            #print("endRowOriginal:" + str(endRowOriginal) + " endColOriginal:" + str(endColOriginal))
            #print("endRowPatch:" + str(endRowPatch) + " endColPatch:" + str(endColPatch))

            # fill the matrix of zeros in the correct position with the tile
            if endRowOriginal > hOrig:
                endRowPatch = endRowPatch - (endRowOriginal - hOrig)
                endRowOriginal = hOrig
                #print("endRowOriginal:" + str(endRowOriginal) + " endColOriginal:" + str(endColOriginal))
            
            if endColOriginal > wOrig:
                endColPatch = endColPatch - (endColOriginal - wOrig)
                endColOriginal = wOrig

            cropPatch = img[0:(endRowPatch), 0:(endColPatch)]
            #print("cropPatch shape: " + str(np.shape(cropPatch)))
            mergedImg[rOriginal : (endRowOriginal), cOriginal:(endColOriginal)] = cropPatch

    #save the completed slice with the output name and imgNumber .png
    os.chdir(folderOutputMerged)
    name = os.getcwd() + filesep + sampleName + "_slice_" + str(sliceNumber).zfill(4) + fileformat

    #Cellpose returns the code of the cell in the pixel 0=background
    if nbits == 8:
        mergedImg = np.uint8((mergedImg>0)*1)
    else:
        mergedImg = np.uint16((mergedImg>0)*1)

    cv2.imwrite(name, mergedImg)
