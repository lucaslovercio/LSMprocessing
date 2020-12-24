close all; clear all; clc;

%PARAMETERS
source_folder = '';%Slices restored
dest_file = '';%TIFF files containing the 3D segmented volume

%Set the resolution in nanometers for anisotropic scans
resX = ;
resY = ;
resZ = ;

%Processing
scale = resX/resZ;

tiffs = dir([source_folder '*.png*']);%List of files (tiles) corresponding to the channel
NumImgs = size(tiffs,1);%Number of tiles

f = waitbar(0,'1','Name','Downsampling and stacking files...',...
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
    
    full_path_source = [source_folder filesep tiffs(i).name];
    slice_original = imread(full_path_source);
    
    %max(imgItMedianClosed2(:))
    img_downsample = imresize(slice_original,scale,'nearest');
    img_downsample = uint8(img_downsample);
    
    
    if i==1
        imwrite(img_downsample,dest_file, 'Compression', 'LZW');
    else
        imwrite(img_downsample,dest_file, 'Compression', 'LZW', 'writemode', 'append');
    end
    
     %clear slice_original; clear vol_original; clear img_downsample; %Release memory for efficiency
    
end

delete(f);

if canceled
    disp('Cancelled by user!!!');
else
    disp('Finished');
end
