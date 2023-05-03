close all; clear all; clc;

folder = '';
sample = '';


%big samples
valueNeural = 100;
valueMesen = 50;

%Do this after deleting noise in 3D Slicer in Neural
%PARAMETERS
file_volume_Neural=             strcat(folder,sample,'/Step_10_',sample,'_Neural_Masked.tiff');
file_volume_Mesen=              strcat(folder,sample,'/Step_10_',sample,'_Mesen_Masked.tiff');
file_volume_Tissues=            strcat(folder,sample,'/Step_10_',sample,'_Tissues_Masked.tiff');

file_volume_Neural_3DSlicer=    strcat(folder,sample,'/Step_10_',sample,'_Neural_Masked_3DSlicer.tiff');

%folderDest = strcat(folder,sample,'/ManualCorrection/');
folderDest = strcat(folder,sample,'/');
mkdir(folderDest);
file_volume_Neural_corrected=   strcat(folderDest,'Step_10_a_',sample,'_Neural.tiff');
file_volume_Mesen_corrected=    strcat(folderDest,'Step_10_a_',sample,'_Mesen.tiff');
file_volume_Tissues_corrected=  strcat(folderDest,'Step_10_a_',sample,'_Tissues.tiff');

%-----------------------------------------------------

%Correct Neural, which is delete in Neural
[volume_Neural,~] = functionReadTIFFMultipage(file_volume_Neural);
[volume_Neural_3DSlicer,~] = functionReadTIFFMultipage(file_volume_Neural_3DSlicer);

volume_Neural_corrected = zeros(size(volume_Neural));
volume_Neural_corrected(volume_Neural_3DSlicer>0) = volume_Neural(volume_Neural_3DSlicer>0);

functionSaveTIFFMultipage(volume_Neural_corrected,file_volume_Neural_corrected,8);

%Correct Neural, which is change label from Neural to Mesen

diffVoxels = xor(volume_Neural>0,volume_Neural_corrected>0);

[volume_Mesen_corrected,~] = functionReadTIFFMultipage(file_volume_Mesen);

%extract label value
labelValue = max(volume_Mesen_corrected(:));
volume_Mesen_corrected(diffVoxels) = labelValue;
functionSaveTIFFMultipage(volume_Mesen_corrected,file_volume_Mesen_corrected,8);

%modify the tissues volume
[volume_Tissues_corrected,~] = functionReadTIFFMultipage(file_volume_Tissues);
labelValue = valueMesen;
volume_Tissues_corrected(diffVoxels) = labelValue;%change the value in those voxels
functionSaveTIFFMultipage(volume_Tissues_corrected,file_volume_Tissues_corrected,8);
