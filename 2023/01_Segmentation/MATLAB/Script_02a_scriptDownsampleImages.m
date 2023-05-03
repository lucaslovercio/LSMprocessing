clear all; close all; clc;
folder = '';
folderOrig = '';
folderDest = '';
ending = '.png';


listing = dir([folderOrig '/' '*.png']);
nFiles = size(listing,1)
imgSize = 256;

showMasks = false;

for i=1:nFiles
    fullpathOrig = strcat(listing(i).folder,filesep,listing(i).name);
    imgRGB = imread(fullpathOrig); sizeRGB = size(imgRGB);
    img = imgRGB(:,:,1);
    img = imresize(img,[imgSize imgSize],'nearest');
    
    %lenghtOrig = strlength(listing(i).name);
    %lenghtEnding = strlength(ending);
    %firstPart = substr(listing(i).name,1,lenghtOrig-lenghtEnding)
    firstPart = listing(i).name;
    
    imwrite(img,strcat(folderDest,filesep,firstPart));
    
end
