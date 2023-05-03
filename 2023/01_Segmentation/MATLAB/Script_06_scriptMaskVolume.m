close all; clear all; clc;

%PARAMETERS
file_volume_Mask=   '...tiff';
file_volume_ToMask= '...tiff';
file_dest_Masked=   '...tiff';

%-----------------------------------------------------

[volumeMask,fileInfoMask] = functionReadTIFFMultipage(file_volume_Mask);
[volumeToMask,fileInfoToMask] = functionReadTIFFMultipage(file_volume_ToMask);

volumeMasked = zeros(size(volumeToMask));
volumeMasked(volumeMask>0) = volumeToMask(volumeMask>0);

functionSaveTIFFMultipage(volumeMasked,file_dest_Masked,8);
