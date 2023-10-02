clear all; close all; clc;
folderVolume =      '';
folderDestTiles =   '';
sampleName =        '';
outputName = strcat(folderDestTiles,filesep,sampleName);

%output parameters
%nbits = 8; %can be 16
nChannels = 1;%can be 3
fileformat = '.png';

patchSize = 1024;%patchSize = 512;%patchSize = 1024;
spx = 100; %neighbourhood
showTiles = false;

%------------------------------------------------------------

listing = dir(folderVolume);
nFiles = length(listing);

usefulSize = patchSize - (2 * spx);


for imgNumber=3:nFiles
    disp(strcat('Processing:',num2str(imgNumber)))
    
    %DEBUG CORRECTION!! THE EXPORT OF TIFFS DO NOT ADD ZERO-PADDING SO THE
    %ORDER ABOVE 1000 SLICES IS WRONG!
    zNumberStr = extractBetween(listing(imgNumber).name,"_Z","_C");
    sliceNumber = str2num(zNumberStr{1});
    %disp(zNumber)
    
    fullpathOrig = strcat(listing(imgNumber).folder,filesep,listing(imgNumber).name);
    img = imread(fullpathOrig); img = img(:,:,1);
    [hOrig,wOrig] = size(img);
    
    imgNorm = functionPercNorm(img);
    
    hImgAugmented = hOrig + patchSize;
    wImgAugmented = wOrig + patchSize;
    imgAugmented =  zeros(hImgAugmented,wImgAugmented);
    imgAugmented((spx + 1):(spx + hOrig),(spx + 1):(spx + wOrig)) = imgNorm;
    
    %figure('Name','Original image');imshow(imgNorm)
    %figure('Name','Augmented image');imshow(imgAugmented)
    
    %Saving the tiles
    ca = {};
    plotIndex = 1;
    numBlocksR = 0;
    for r = 1 : usefulSize : hOrig %iteration inside the true image
        numBlocksC = 0;
        numBlocksR = numBlocksR + 1;
        for c = 1: usefulSize : wOrig
            numBlocksC = numBlocksC + 1;
            %Add neighbourhood to tile
            
            %top left corner
            r1 = r;
            c1 = c;
            
            %Upper neighbourhood in the true image
            r0 = r1 - spx;
            c0 = c1 - spx;
            
            %bottom right neighbourhood in the true image
            r2 = r0 + patchSize - 1;
            c2 = c0 + patchSize - 1;
            
            %traslation to augmented image
            r0_a = r0 + spx;
            c0_a = c0 + spx;
            r2_a = r2 + spx;
            c2_a = c2 + spx;
            
            %disp(strcat('r1',num2str(r1),'-c1',num2str(c1),'-r0',num2str(r0),'-c0',num2str(c0),...
            %    'r0_a',num2str(r0_a),'-c0_a',num2str(c0_a),'-r2_a',num2str(r2_a),'-c2_a',num2str(c2_a)));
            
            block = imgAugmented(r0_a:r2_a,c0_a:c2_a);
            
            name = strcat(outputName,'_slice_',num2str(sliceNumber,'%04.f'),...
                '_block_',num2str(numBlocksR,'%04.f'),'_',num2str(numBlocksC,'%04.f'),'.png');
            if nChannels==3
                rgb = cat(3,block,block,block);
            else
                rgb = block;
            end
            imwrite(rgb,name);
            
            plotIndex = plotIndex + 1;
            ca{numBlocksR,numBlocksC} = block;
            
        end
        
    end
    
    
    if showTiles
        hFig = figure('Name',listing(imgNumber).name);
        plotIndex = 1;
        for r = 1 : numBlocksR
            for c = 1 : numBlocksC
                %fprintf('plotindex = %d,   c=%d, r=%d\n', plotIndex, c, r);
                % Specify the location for display of the image.
                figure(hFig);subplot(numBlocksR, numBlocksC, plotIndex);
                % Extract the numerical array out of the cell
                % just for tutorial purposes.
                rgbBlock = ca{r,c};
                imshow(rgbBlock); % Could call imshow(ca{r,c}) if you wanted to.
                [rowsB, columnsB, numberOfColorBandsB] = size(rgbBlock);
                % Make the caption the block number.
                caption = sprintf('Block #%d of %d\n%d rows by %d columns', ...
                    plotIndex, numBlocksR*numBlocksC, rowsB, columnsB);
                title(caption);
                drawnow;
                % Increment the subplot to the next location.
                plotIndex = plotIndex + 1;
            end
        end
    end
    
    %Save slice data
    infoName = strcat(outputName,'_slice_',num2str(sliceNumber,'%04.f'),'_info');
    save(strcat(infoName,'.mat'),...
        'sampleName','sliceNumber','numBlocksR','numBlocksC','imgNumber',...
        'hOrig','wOrig','fullpathOrig','patchSize','hImgAugmented', 'wImgAugmented',...
        'usefulSize','spx');
    
    fid = fopen(strcat(infoName,'.txt'),'wt');
    fprintf(fid, strcat(fullpathOrig,'\n'));
    fprintf(fid, strcat('sampleName:',sampleName,'\n'));
    fprintf(fid, strcat('imgNumber:',num2str(sliceNumber),'\n'));
    fprintf(fid, strcat('hOrig:',num2str(hOrig),'\n'));
    fprintf(fid, strcat('wOrig:',num2str(wOrig),'\n'));
    fprintf(fid, strcat('hImgAugmented:',num2str(hImgAugmented),'\n'));
    fprintf(fid, strcat('wImgAugmented:',num2str(wImgAugmented),'\n'));
    fprintf(fid, strcat('patchSize:',num2str(patchSize),'\n'));
    fprintf(fid, strcat('numBlocksR:',num2str(numBlocksR),'\n'));
    fprintf(fid, strcat('numBlocksC:',num2str(numBlocksC),'\n'));
    fprintf(fid, strcat('spx:',num2str(spx),'\n'));
    fprintf(fid, strcat('usefulSize:',num2str(usefulSize),'\n'));
    
    fclose(fid);
    
end
