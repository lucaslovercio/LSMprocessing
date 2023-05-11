close all; clear all; clc;

previous = 'E11.0';
%next = 'E11.5';

%nexts_names = {'Dec2_E10_9','Dec2_E10_9_flipped','Dec2_E10_10','Dec2_E10_10_flipped','Dec2_E10_17','Dec2_E10_17_flipped',...
%     'Dec2_E105_3','Dec2_E105_3_flipped','Oct27_E105_4', 'Oct27_E105_4_flipped'};
%nexts_names = {'Feb3_E11_4_v2','Feb3_E11_4_v2_flipped','Feb12_E115_1','Feb12_E115_1_flipped',...
% 'Feb12_E115_2','Feb12_E115_2_flipped','Feb12_E115_8','Feb12_E115_8_flipped','May10_E11_2','May10_E11_2_flipped'};
nexts_names = {'Feb12_E115_3','Feb12_E115_3_flipped','Feb12_E115_4','Feb12_E115_4_flipped','Feb12_E115_5','Feb12_E115_5_flipped',...
     'Feb12_E115_7','Feb12_E115_7_flipped','May10_E11_8','May10_E11_8_flipped'};

only_nonzero = false;
only_shapes_involved = true;
 
nFiles = length(nexts_names);
vCorr = double(zeros(nFiles,1));
vPropShapeChange = double(zeros(nFiles,1));
for j=1:nFiles
    next = nexts_names{j};
    disp('------------------------');
    disp(next);
    folder =        ''; %Folder with all samples registered to the same objective
    pathVolumeDifferenceNifti = strcat(folder,'Diff_volumes_not_masked_noW',filesep,...
        'VolumeDifference_',previous,'_Tissues_to_',next,'_mesen.nii');
    pathProliferation = strcat(folder,previous,filesep,previous,'_PHH3_Smooth_w_30_HistNorm.tiff');
    pathMask = strcat(folder,'MesenMajority03_cropped_shifted_120_cropped.tiff');
    
    %Read volumes
    diff_volume = niftiread(pathVolumeDifferenceNifti);
    prolif_volume = functionReadTIFFMultipage(pathProliferation);
    mask_volume  = functionReadTIFFMultipage(pathMask); mask_volume = mask_volume>0;
    
    if only_nonzero
        non_zero = abs(diff_volume)>0;
        disp('nVoxels mask and non-zero diff voxels');
        sum_mask = sum(mask_volume(:)); disp(sum_mask);
        disp(sum(non_zero(:)));
        mask_volume = mask_volume & non_zero;
        disp('mask and non-zero');
        sum_mask_nonzero = sum(mask_volume(:)); disp(sum_mask_nonzero);
        disp('Prop shape change');
        vPropShapeChange(j)=sum_mask_nonzero/sum_mask;
        clear sum_mask
        disp(vPropShapeChange(j));
    end
    
    if only_shapes_involved
        
    
        mean_shape_prev_path = strcat(folder,previous,filesep,previous,'_Tissues_Affine_majority.tiff');
        mean_shape_next_path =	strcat(folder,next,'_toE10.5/',next,'_Tissues_Masked_toE10.5.tiff');
        mean_shape_prev = functionReadTIFFMultipage(mean_shape_prev_path);
        mean_shape_next = functionReadTIFFMultipage(mean_shape_next_path);
        mask_volume = mask_volume & ((mean_shape_prev>0) |  (mean_shape_next>0));
        clear mean_shape_prev mean_shape_next
    end
    
    %Values in the mask
    diff_values = double(diff_volume(mask_volume));
    prolif_values = double(prolif_volume(mask_volume));
    clear mask_volume diff_volume prolif_volume
    
    %Correlation
    matCorr = horzcat(prolif_values,diff_values);
    CorrMatrix = corr(matCorr);
    disp(next);
    disp(CorrMatrix(2,1));
    vCorr(j)=CorrMatrix(2,1);
end




%Mean shapes, Including zeros (no volumen change)
%E10.0-E10.5 : -0.5389
%E10.5-E11.0 : -0.2982
%E11.0-E11.5 : 0.0016

%Excluding no volume change (<5000 vol diff)
%E10.0-E10.5 : -0.4411
%E10.5-E11.0 : -0.2372
%E11.0-E11.5 : 0.0025

%-------------------------------------------------------------------
%Including zeros (no volumen change)
%E10.0-E10.5 :
%Dec2_E10_9: -0.2584
%Dec2_E10_9_flipped: -0.2655
%Dec2_E10_10: -0.3820
%Dec2_E10_10_flipped: -0.3861
%Dec2_E10_17: -0.3981
%Dec2_E10_17_flipped: -0.3880
%Dec2_E105_3: -0.4187
%Dec2_E105_3_flipped: -0.4003
%Oct27_E105_4: -0.4469
%Oct27_E105_4_flipped: -0.4455
 
%E10.5-E11.0 :
%Feb3_E11_4_v2: -0.2449
%Feb3_E11_4_v2_flipped: -0.2541 
%Feb12_E115_1: -0.1698
%Feb12_E115_1_flipped: -0.1741 
%Feb12_E115_2: -0.2158
%Feb12_E115_2_flipped: -0.2036 
%Feb12_E115_8: -0.2644
%Feb12_E115_8_flipped: -0.2677 
%May10_E11_2: -0.2460
%May10_E11_2_flipped: -0.2487 

%E11.0-E11.5 :
%Feb12_E115_3: 0.0090
%Feb12_E115_3_flipped: 0.0043
%Feb12_E115_4: 0.0265
%Feb12_E115_4_flipped: 0.0195
%Feb12_E115_5: 0.0196
%Feb12_E115_5_flipped: 0.0133
%Feb12_E115_7: 0.0094
%Feb12_E115_7_flipped: 0.0050
%May10_E11_8: -0.0216
%May10_E11_8_flipped: -0.0200

%-------------------------------------------------------------------
%Excluding no volume change (<800 vol diff)
%E10.0-E10.5 :
%Dec2_E10_9: -0.3037
%Dec2_E10_9_flipped: -0.3076
%Dec2_E10_10: -0.4048
%Dec2_E10_10_flipped: -0.4137
%Dec2_E10_17: -0.4387
%Dec2_E10_17_flipped: -0.4272
%Dec2_E105_3: -0.4378
%Dec2_E105_3_flipped: -0.4232
%Oct27_E105_4: -0.4538
%Oct27_E105_4_flipped: -0.4541

%E10.5-E11.0 :
%Feb3_E11_4_v2: 
%Feb3_E11_4_v2_flipped: 
%Feb12_E115_1: 
%Feb12_E115_1_flipped: 
%Feb12_E115_2: 
%Feb12_E115_2_flipped: 
%Feb12_E115_8: 
%Feb12_E115_8_flipped: 
%May10_E11_2: 
%May10_E11_2_flipped: 

%-------------------------------------------------------------------
%Excluding no volume change (<5000 vol diff)
%E10.0-E10.5 : 
%Dec2_E10_9: -0.4095
%Dec2_E10_9_flipped: -0.4102
%Dec2_E10_10: -0.4480
%Dec2_E10_10_flipped: -0.4697
%Dec2_E10_17: -0.4702
%Dec2_E10_17_flipped: -0.4707
%Dec2_E105_3: -0.4731
%Dec2_E105_3_flipped: -0.4664
%Oct27_E105_4: -0.4709
%Oct27_E105_4_flipped: -0.4756

%E10.5-E11.0 :
%Feb3_E11_4_v2: -0.2703
%Feb3_E11_4_v2_flipped: -0.2665 
%Feb12_E115_1: -0.2807
%Feb12_E115_1_flipped: -0.2864 
%Feb12_E115_2: -0.3438
%Feb12_E115_2_flipped: -0.3154 
%Feb12_E115_8: -0.2823
%Feb12_E115_8_flipped: -0.2692 
%May10_E11_2: -0.2546
%May10_E11_2_flipped: -0.2531 

%E11.0-E11.5 :
%Feb12_E115_3: 0.0285
%Feb12_E115_3_flipped: 0.0240
%Feb12_E115_4: -0.0320
%Feb12_E115_4_flipped: -0.0295
%Feb12_E115_5: 0.0192
%Feb12_E115_5_flipped: 0.0134
%Feb12_E115_7: 0.0028
%Feb12_E115_7_flipped: 0.0170
%May10_E11_8: -0.0114
%May10_E11_8_flipped: -0.0134
