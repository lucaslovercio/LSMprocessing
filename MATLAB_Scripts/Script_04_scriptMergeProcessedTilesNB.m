clear all; close all; clc;
folderTilesOriginal = '';%Folder with unprocessed tiles, to recover the information from .mat files
folderTilesProcessed = '';%Folder with segmented tiles using any U-net
folderOutputMerged = '';%Folder where the restored segmented slices will be saved


%output parameters
nbits = 8; %can be 16
fileformat = '.png';

%------------------------------------------------

%Look for all the .mat where information is stored
filesToSearch = [folderTilesOriginal filesep '*.mat'];
listing = dir(filesToSearch);
nFilesMat = length(listing)

% for each file in the mat lis
for i=1:nFilesMat
    
    fullpathMat = strcat(listing(i).folder,filesep,listing(i).name);
    load(fullpathMat);
    %read the output name and imgNumber
    
    %build the first part of the string
    
    firstPatPath = strcat(folderTilesProcessed,filesep,sampleName,'_slice_',num2str(sliceNumber,'%04.f'));
    disp(strcat('Processing:',num2str(sliceNumber)))
    %read the original hOrig y wOrig, create an empty matrix of zeros
    if nbits == 8
        mergedImg = uint8(zeros(hOrig,wOrig));
    else
        mergedImg = uint16(zeros(hOrig,wOrig));
    end
    %read the numBlocksR and numBlocksC
    for r = 1 : numBlocksR
        for c = 1 : numBlocksC
            %for r in Rows and C in Cols
            %reconstruct the png filename
            processedFile = strcat(firstPatPath,'_block_',num2str(r,'%04.f'),'_',num2str(c,'%04.f'),'.png');
            %read the file
            imgRGB = imread(processedFile); img = imgRGB(:,:,1);
            
            [hProcessed,wProcessed] = size(img);
            if hProcessed~= patchSize
                img = imresize(img,[patchSize,patchSize],'nearest');
            end
            
            img = img(spx + 1:end-spx,spx + 1:end-spx);
            
            infoProcessed = imfinfo(processedFile);
            if r == 1 && c==1 %first iteration
                if infoProcessed.BitDepth == 16
                    mergedImg = uint16(mergedImg);
                else
                    mergedImg = uint8(mergedImg);%Only one 8bits channel!!!
                end
            end
                
            %multiply the r and c per the patchSize
            r2 = (r-1) * usefulSize + 1;
            c2 = (c-1) * usefulSize + 1;
            
            endRow = r2 + usefulSize -1;
            endCol = c2 + usefulSize -1;
            endRowPatch = usefulSize;
            endColPatch = usefulSize;
            
            %fill the matrix of zeros in the correct position with the tile
            if endRow > hOrig
                endRowPatch = usefulSize - (endRow - hOrig);
                endRow = hOrig;
            end
            if endCol > wOrig
                endColPatch = usefulSize - (endCol - wOrig);
                endCol = wOrig;
            end
            
            mergedImg(r2:endRow,c2:endCol) = img(1:endRowPatch,1:endColPatch);
            
        end
    end
    
    %save the completed slice with the output name and imgNumber .png
    imwrite(mergedImg,strcat(folderOutputMerged,filesep,sampleName,'_slice_',num2str(sliceNumber,'%04.f'),fileformat));
    
end
