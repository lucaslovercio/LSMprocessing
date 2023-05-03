close all; clear all; clc;

%PARAMETERS
source_folder = '';
dest_file =     '....tiff';

%In nanometers!!!
%Search the voxel dimensions in Superhero Embryos!!!!, Sheet 2
resX = 913.89;
resY = 913.89;
resZ = 4940;

%Processing
scale = resX/resZ;

tiffs = dir([source_folder '*.png*']);%List of files (tiles) corresponding to the channel
NumImgs = size(tiffs,1)%Number of tiles

for i=1:NumImgs
    
    disp(strcat('Processing:',num2str(i)))
    
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
