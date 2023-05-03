close all; clear all; clc;

%PARAMETERS
source_folder_segmented =   '';
dest_file_Neural=           '...tiff';
dest_file_Mesen=            '...tiff';

%1024 side
valueNeural = 50;
valueMesen = 100;

%In nanometers!!!

resX = 913.89;
resY = 913.89;
resZ = 4940;


%------------------------------------------------------------------

%Processing
scale = resX/resZ;

tiffs = dir([source_folder_segmented '*.png*']);%List of files (tiles) corresponding to the channel
NumImgs = size(tiffs,1);%Number of tiles

f = waitbar(0,'1','Name','Filtering DAPI, downsampling and stacking files...',...
    'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
canceled = false;

for i=1:NumImgs
    
    if getappdata(f,'canceling')
        canceled = true;
        delete(f);
        break
    end
    % Update waitbar and message
    waitbar(i/NumImgs,f,sprintf('%i / %i',i,NumImgs))
    
    full_path_source = [source_folder_segmented filesep tiffs(i).name];
    slice_original = imread(full_path_source);
    
    %tissue thresholding
    sliceNeural = slice_original == valueNeural;
    sliceMesen = slice_original == valueMesen;
    
    img_downsample_mesen = imresize(sliceMesen,scale,'nearest');
    img_downsample_mesen = uint8(img_downsample_mesen);
    
    img_downsample_neural = imresize(sliceNeural,scale,'nearest');
    img_downsample_neural = uint8(img_downsample_neural);
    
    if i==1
        imwrite(img_downsample_neural,dest_file_Neural, 'Compression', 'LZW');
    else
        imwrite(img_downsample_neural,dest_file_Neural, 'Compression', 'LZW', 'writemode', 'append');
    end
    
    if i==1
        imwrite(img_downsample_mesen,dest_file_Mesen, 'Compression', 'LZW');
    else
        imwrite(img_downsample_mesen,dest_file_Mesen, 'Compression', 'LZW', 'writemode', 'append');
    end
    
     %clear slice_original; clear vol_original; clear img_downsample; %Release memory for efficiency
    
end

delete(f);

if canceled
    disp('Cancelled by user!!!');
else
    disp('Finished');
end
