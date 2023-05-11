clear all; close all; clc;

rootdir =      '';
pathAll =      '...tiff'; %TIFF stack to generate

filelist = dir(fullfile(rootdir, '**/*....tiff')); %Volumes to sum

nVolumes = length(filelist)

volumeTemp = functionReadTIFFMultipage(strcat(filelist(1).folder,filesep,filelist(1).name));
volumeSum = double(zeros(size(volumeTemp)));
clear volumeTemp
for i=1:nVolumes
    path1 = strcat(filelist(i).folder,filesep,filelist(i).name);
    disp(filelist(i).name);
    
    volume1 = functionReadTIFFMultipage(path1);
    
    volumeSum = volumeSum + double(volume1);
    
    clear volume1
     
end

functionSaveTIFFMultipage(volumeSum./double(nVolumes),pathAll,8);%As the volumes are normalized 0-255, this is the avg