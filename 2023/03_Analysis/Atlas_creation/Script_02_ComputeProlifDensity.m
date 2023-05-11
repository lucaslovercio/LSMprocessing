close all; clear all; clc;

%This convolution can be applied to obtain the proliferation map of one
%sample of an atlas if the proliferating volume is the sum of all of the
%group

%PARAMETERS


age = 'E11.5';
histNorm = true;
regionSizes = [30];

for regionSize = regionSizes

    
    rootdir = ''; %Folder where to search recursively
    ending = '_PHH3_sum.tiff'; % Proliferating cells volumes ending or total sum of proliferation in a group
    filelist = dir(fullfile(rootdir, strcat('**/*',ending)));
    
    path_volume_shape = ''; %Shape atlas
   
    nVolumes = length(filelist)
    %nVolumes = 3;

    for i=1:nVolumes
        disp(filelist(i).name);
        sampleName1 = split(filelist(i).name,ending); sampleName = sampleName1{1};
        disp(sampleName)
        path_volume_PHH_segmentation =  strcat(filelist(i).folder,filesep,filelist(i).name);


        dest_file_PHHSmoothed =     strcat(filelist(i).folder,filesep,sampleName,'_PHH3_Smooth','_w_',num2str(regionSize),'.tiff');
        dest_file_PHHSmoothedNorm = strcat(filelist(i).folder,filesep,sampleName,'_PHH3_Smooth_Norm','_w_',num2str(regionSize),'.tiff');

        if histNorm
            dest_file_PHHSmoothed =     strcat(filelist(i).folder,filesep,sampleName,'_PHH3_Smooth','_w_',num2str(regionSize),'_HistNorm.tiff');
            dest_file_PHHSmoothedNorm = strcat(filelist(i).folder,filesep,sampleName,'_PHH3_Smooth_Norm','_w_',num2str(regionSize),'_HistNorm.tiff');
        end

        %-----------------------------------------------------------------------------------
        %Script

        infoTIFFCellPHH  = imfinfo(path_volume_PHH_segmentation);
        nSlices = length(infoTIFFCellPHH);
        slice = imread(path_volume_PHH_segmentation, 1, 'Info', infoTIFFCellPHH);
        [h,w] = size(slice);

        volumePHH = uint32(zeros(h,w,nSlices));%Change the data type to the proper path_volume_cell_segmentation data type

        infoTIFFShape  = imfinfo(path_volume_shape);
        volume_shape = uint8(zeros(h,w,nSlices));

        for i=1:nSlices

            sliceShape = uint8(imread(path_volume_shape, i, 'Info', infoTIFFShape));
            volume_shape(:,:,i) = sliceShape;

            slicePHH = uint32(imread(path_volume_PHH_segmentation, i, 'Info', infoTIFFCellPHH));
            volumePHH(:,:,i) = slicePHH;

        end

        B = ones(regionSize,regionSize,regionSize);
        disp('Convoluting....');

        volumePHHSmoothed = convn(volumePHH,B,'same');
        disp('Convolution Finished');


        volumePHHSmoothed_corrected = zeros(size(volumePHHSmoothed));
        volumePHHSmoothed_corrected(volume_shape>0) = volumePHHSmoothed(volume_shape>0);

        if histNorm
            volumePHHSmoothed_corrected_Norm = functionWeigertNorm(volumePHHSmoothed_corrected);
        else
            volumePHHSmoothed_corrected_Norm = double(volumePHHSmoothed_corrected)./ double(max(volumePHHSmoothed_corrected(:)));
        end


        for i=1:nSlices

            sliceToSavePHH = volumePHHSmoothed_corrected(:,:,i);
            sliceToSavePHHNorm = volumePHHSmoothed_corrected_Norm(:,:,i);
            if i==1

                imwrite(sliceToSavePHH,         dest_file_PHHSmoothed,      'Compression', 'LZW');
                imwrite(sliceToSavePHHNorm,     dest_file_PHHSmoothedNorm,  'Compression', 'LZW');

            else

                imwrite(sliceToSavePHH,         dest_file_PHHSmoothed,      'Compression', 'LZW', 'writemode', 'append');
                imwrite(sliceToSavePHHNorm,     dest_file_PHHSmoothedNorm,  'Compression', 'LZW', 'writemode', 'append');

            end

            %clear slice_original; clear vol_original; clear img_downsample; %Release memory for efficiency

        end
        
        %Showing new histogram for proliferation
        position = [50 50 800 800];

        vForHistogram = volumePHHSmoothed_corrected(:);
        vForHistogramNonzero = vForHistogram(vForHistogram>0);
        fHistNonZero = figure('Name','Histogram PHH NonZero','Position',position);
        histogram(vForHistogramNonzero);title(strrep(strcat(dest_file_PHHSmoothed,'_histAfterCorrection.png'),'_','-'));
        saveas(fHistNonZero,strcat(dest_file_PHHSmoothed,'_histAfterCorrection.png'));

        close all;
    end
end