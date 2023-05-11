close all; clear all; clc;

%PARAMETERS
age = 'E11.0';
shiftValue = 120; %Depth of mask, as different were tested

rootdir = '';
front_mask_path = strcat(rootdir,age,'/LM38RigidReg_itk/MeanNNTransform/',age,'_erode_and_mask_',num2str(shiftValue),'.tiff'); %Frontal mask location

folderAffine = strcat(rootdir,age,'/LM38RigidReg_itk/MeanNNTransform'); %Folder with affine registered volumes
ending = '_PHH3_Smooth_Norm_w_30_HistNorm.tiff'; % Ending of proliferation maps
endingMat   = strcat(ending,'_proliferationValuesFront_',num2str(shiftValue),'.mat'); %Saving the proliferation in the mask as mat for PLS
endingTIFF  = strcat(ending,'_proliferationValuesFront_',num2str(shiftValue),'.tiff'); %Saving the proliferation in the mask as mat for visualization
filelistProliferation = dir(fullfile(folderAffine, strcat('**/*',ending)));

nVolumes = length(filelistProliferation)

maskFront = functionReadTIFFMultipage(front_mask_path);
maskFront = maskFront>0;%To make if a binary mask


for i=1:nVolumes
    disp('-------')
    disp(filelistProliferation(i).name);
    path_volume_Proliferation =  strcat(filelistProliferation(i).folder,filesep,filelistProliferation(i).name);
    
    sampleName = split(filelistProliferation(i).name,ending); sampleName = sampleName{1};

    proliferationNorm = functionReadTIFFMultipage(path_volume_Proliferation);
    maxProlif = max(proliferationNorm(:)); disp(strcat('max prolif: ',num2str(maxProlif))); %should be 255
    minProlif = min(proliferationNorm(:)); disp(strcat('min prolif: ',num2str(minProlif))); %should be 0
   
    
    proliferationNormDouble = double(proliferationNorm)./double(maxProlif);
        
    valuesProliferation = proliferationNormDouble(maskFront);
    disp(num2str(mean(valuesProliferation)));
    matPath = strcat(filelistProliferation(i).folder,filesep,sampleName,endingMat);
    save(matPath,'valuesProliferation','sampleName');
    
    clear valuesProliferation proliferationNormDouble;
    
    proliferationInMaskVolume = uint8(zeros(size(proliferationNorm)));
    proliferationInMaskVolume(maskFront) = proliferationNorm(maskFront);
    volumeOutputPath = strcat(filelistProliferation(i).folder,filesep,sampleName,endingTIFF);
    disp(volumeOutputPath);
    functionSaveTIFFMultipage(proliferationInMaskVolume,volumeOutputPath,8);
    
end

