%% Variable Needed
% Calibration data
% sfX : Image shift along X
% sfY : Image shift along Y
% imgRaw : HDR cannula image
clear; clc; %close all;
%load get calibibration file
%path = 'G:\plant\Red_stain\plant_crop\';

% load([path,'b.mat']);
% background = frames.imgRaw;

load(['20230821_PSF_noNormalization.mat']);
%load(['20230731.mat']);
%%
%load test image CCM raw
%load(['/media/alingold/MenonLab/20230807_50um_25cal_100umWide/x24_y24_z0.mat'],'objCCM', 'objRef');
%load(["C:\Cannula Microscope\20230822_5x5_multiBeadTest\x4_y0_z0.mat"],'objCCM', 'objRef');
%load(["C:\Cannula Microscope\20230822_plantTest\9.mat"],'objCCM', 'objRef');
load(["C:\Cannula Microscope\20230823_plantTest\x0_y0_z0.mat"],'objCCM', 'objRef');
%load(['/media/alingold/MenonLab/20230803_50um_bead_calTest/x1_y0_z0.mat'],'objCCM', 'objRef');
%load(['/media/alingold/MenonLab/20230731_50um_25cal/x20_y10_z0.mat'],'objCCM', 'objRef');
%load(['/media/alingold/MenonLab/CCM Ruipeng/2022_02_21_psf_for_plant_1/x15_y15_z0.mat'],'imgHDR', 'objRef');
% 20230713 REFcrop = [400, 840, 440, 860];

REFcrop = [420, 1050, 220, 900]
%scalefactor = 0.4;
%image = imresize(imgHDR, scalefactor);
%image = imgHDR;
% apply mask
%load(['Mask.mat'], 'mask')
size(mask)
size(objCCM)
% % Finding the min and max value within the mask
% min_value_within_mask = min(image(mask));
% max_value_within_mask = max(image(mask));
% 
% % Normalizing the image
% image = (image - min_value_within_mask) / (max_value_within_mask - min_value_within_mask);
imgHDR = objCCM;
imgHDR(313, 457) = 0;
masked_values = imgHDR .* mask;
min_value_within_mask = min(masked_values(masked_values ~= 0))
max_value_within_mask = max(max(masked_values))
A_reshaped = reshape(masked_values, [], 1);
sorted_values = sort(A_reshaped, 'descend');
max_value_within_mask = sorted_values(2);
% Normalizing the image
imgHDR = (imgHDR - min_value_within_mask) / (max_value_within_mask - min_value_within_mask);

% Applying the mask
masked_image = imgHDR .* mask;
% crop to size
cropped_image = masked_image(top_left_row:bottom_right_row, top_left_col:bottom_right_col);
% resize to 128x128
image = imresize(cropped_image, [imageSize imageSize]);


tic;
for i = 1:10

Ref = objRef; % pull ref img (ground truth)
%Recon = frames.objRaw;
lambda = 10000;
frames.imgPro = reshape(image, [1, imageSize^2])';

% if (frames.sfX ~= 0 || frames.sfY ~= 0) % this is for shifting, not used
%     refPro = imref2d(size(frames(f).imgPro));
%     tform  = affine2d([1 0 0; 0 1 0; frames.sfX frames.sfY 1]);
%     frames.imgPro = imwarp(frames.imgPro,tform,'OutputView',refPro,'FillValues',0);
% end
% frames.imgPro = imresize(frames.imgPro,size(imgFOV)); % resize to FOV

% if size(frames.imgPro) ~= size(imgFOV)
%     error('Image size does not match with calibration. Check your image subregion setting.');
% end
% apply process to example image which is save in calibration file
frames.imgPro = single(frames.imgPro) - DC;
minImg = min(frames.imgPro(imgFOV));
%minImg = min(frames.imgPro);
frames.imgPro(~imgFOV) = minImg;
%frames.imgPro = minImg;
frames.imgPro = double(frames.imgPro);

clear refPro tform minImg;

%y = double(frames.imgPro(imgFOV));
y = double(frames.imgPro);
if (min(lambda)<0)
  error('Illegal regularization parameter lambda')
end
% tikhonov     
n = size(V,1);    
beta = U'*y;
zeta = s(:,1).*beta;
ll = length(lambda);
x_lambda = zeros(n,ll);
    
% Treat each lambda separately.
for i=1:ll
    x_lambda(:,i) = V*(zeta./(s.^2 + lambda(i)^2));
end
% extract guess from x_lambda
%objRaw = zeros(1,psfSize^2*zSize);
%objFOVInd = find(objFOV==1);
%objRaw(objFOVInd) = x_lambda;
objRaw = x_lambda;
objRaw = reshape(objRaw,[psfSize psfSize zSize]);
objRaw = abs(objRaw);
end
toc

threshold = 0;
for n = 1:size(objRaw)
    for m = 1:size(objRaw)
        if objRaw(n, m) < threshold
            objRaw(n, m) = 0;
        end
    end
end


REFcrop = [420, 1050, 220, 900];
Ref_cropped = Ref(REFcrop(3):REFcrop(4), REFcrop(1):REFcrop(2));
Ref_cropped = rot90(Ref_cropped);
Ref_cropped_scaled = imresize(Ref_cropped, [100, 100]);
REFcrop2 = [35 90 45 100];
Ref_cropped_2 = Ref_cropped_scaled(REFcrop2(3):REFcrop2(4), REFcrop2(1):REFcrop2(2));
objRaw_cropped = objRaw(REFcrop2(3):REFcrop2(4), REFcrop2(1):REFcrop2(2));




%Ref_cropped_2 = rot90(Ref_cropped_2);

figure;

% subplot(2,3,1);
% imagesc(Ref);
% colormap gray;
% title('uncropped ground truth');
% 
% subplot(2,3,2);
% imagesc(Ref_cropped);
% colormap gray;
% title('cropped ground truth to scan');

subplot(2,2,1);
imagesc(Ref_cropped_2);
colormap gray;
title('cropped ground truth to FOV');

subplot(2,2,2);
imagesc(image);
colormap gray;
title('Processed Fiber');

subplot(2,2,3);
imagesc(objRaw);
colormap gray;
title('reconstruction');

subplot(2,2,4);
imagesc(objRaw_cropped);
colormap gray;
title('reconstruction cropped to FOV');
