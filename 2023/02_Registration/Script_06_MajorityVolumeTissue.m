clear all; close all; clc;

rootdir             = '...';
pathMajorityTissue  = 'E11.5_Tissues_Affine_majority.tiff';
pathMesen           = 'E11.5_Mesen_Affine_majority.tiff';
pathNeural          = 'E11.5_Neural_Affine_majority.tiff';

filelist1 = dir(fullfile(rootdir, '**/*_Tissues_Masked_toE10.5.tiff'));

listSamples = {'Feb12_E115_3','Feb12_E115_4','Feb12_E115_5','Feb12_E115_7','May10_E11_8'};

nSamples = length(listSamples)
nVolumes = length(filelist1)
idx_filename = 1;
for i=1:nVolumes
    name_temp = filelist1(i).name;
    for j=1:nSamples
        sample_temp = listSamples{j};
        if contains(name_temp,sample_temp)
            filelist(idx_filename).folder = filelist1(i).folder;
            filelist(idx_filename).name = filelist1(i).name;
            idx_filename = idx_filename + 1;
        end
    end
    
end

nVolumes = length(filelist)

volumeTemp = functionReadTIFFMultipage(strcat(filelist(1).folder,filesep,filelist(1).name));
[h,w,z] = size(volumeTemp);
nLabels = 3;%Mesen, Neural, Background
volumeSumBackground = uint8(zeros(h,w,z));
volumeSumMesen = uint8(zeros(h,w,z));
volumeSumNeural = uint8(zeros(h,w,z));

valueNeural = 100;
valueMesen = 50;
valueBackground = 0;



for i=1:nVolumes
    path1 = strcat(filelist(i).folder,filesep,filelist(i).name);
    disp(filelist(i).name)
    
    volume1 = functionReadTIFFMultipage(path1);
    
    volumeSumNeural(volume1==valueNeural) = volumeSumNeural(volume1==valueNeural) + 1;
    volumeSumMesen(volume1==valueMesen) = volumeSumMesen(volume1==valueMesen) + 1;
    volumeSumBackground(volume1==valueBackground) = volumeSumBackground(volume1==valueBackground) + 1;
     
end

volumeSum = uint8(zeros(h,w,z,3));
volumeSum(:,:,:,1) = volumeSumMesen;%channel 1 mesen
volumeSum(:,:,:,2) = volumeSumNeural;%channel 2 neural
volumeSum(:,:,:,3) = volumeSumBackground;%channel 3 background

clear volumeSumBackground volumeSumNeural volumeSumMesen

[~,argmax] = max(volumeSum,[],4);

size(argmax);
min(argmax(:));
max(argmax(:));

%volumeMajority = volumeSum > TH;
volumeOutput = uint8(zeros(h,w,z));
volumeOutput(argmax == 1) = valueMesen;
volumeOutput(argmax == 2) = valueNeural;
volumeOutput(argmax == 3) = valueBackground;

functionSaveTIFFMultipage(volumeOutput,pathMajorityTissue,8);

functionSaveTIFFMultipage(uint8(argmax == 1),pathMesen,8);
functionSaveTIFFMultipage(uint8(argmax == 2),pathNeural,8);
