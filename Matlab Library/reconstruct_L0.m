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

load(['/home/alingold/matlab/20230821_PSF_normFix.mat']);
%load(['20230731.mat']);
%%
%close all
%load test image CCM raw
%load(['/media/alingold/MenonLab/20230807_50um_25cal_100umWide/x24_y24_z0.mat'],'objCCM', 'objRef');
%load(["/media/alingold/MenonLab/20230822_5x5_multiBeadTest/x4_y0_z0.mat"],'objCCM', 'objRef');
load(["/media/alingold/MenonLab/20230901_plant_Exposure1-3000ms/250ms.mat"],'objCCM', 'objRef');
%load(["C:\Cannula Microscope\20230822_plantTest\9.mat"],'objCCM', 'objRef');
%load(["/media/alingold/MenonLab/20230829_cottonFIber/x0_y100_z0.mat"],'objCCM', 'objRef');
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
min_value_within_mask = min(min(masked_values))
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

tic
% Initialize the parameters
Omega = []; % Index set to keep track of selected columns of A
r = y;      % Initial residual
k = 1000;     % Maximum number of non-zero coefficients or stopping criterion
N = size(A, 2); % Number of columns in A

% Begin OMP Algorithm
for iter = 1:k
    % Step 1: Find the index that has maximum correlation with the residual
    [~, idx] = max(abs(A' * r));
    
    % Step 2: Add the found index to the set Omega
    Omega = [Omega, idx];
    
    % Step 3: Update x by solving the least-squares problem
    x_ls = A(:, Omega) \ y;
    
    % Step 4: Update the residual
    r = y - A(:, Omega) * x_ls;
    
    % Check some stopping criterion if you have one
    % For example: if norm(r) < some_tolerance, break; end
end

% Prepare the final sparse x
x_sparse = zeros(N, 1);
x_sparse(Omega) = x_ls;

% You might still want to reshape this vector according to your original image dimensions
objRaw = reshape(x_sparse, [psfSize, psfSize, zSize]);
toc

threshold = 0;
for n = 1:size(objRaw)
    for m = 1:size(objRaw)
        if objRaw(n, m) < threshold
            objRaw(n, m) = 0;
        end
    end
end

Ref = objRef; % pull ref img (ground truth)
REFcrop = [420, 1050, 220, 900];
Ref_cropped = Ref(REFcrop(3):REFcrop(4), REFcrop(1):REFcrop(2));
Ref_cropped = rot90(Ref_cropped);
Ref_cropped_scaled = imresize(Ref_cropped, [100, 100]);
REFcrop2 = [35 90 45 100];
Ref_cropped_2 = Ref_cropped_scaled(REFcrop2(3):REFcrop2(4), REFcrop2(1):REFcrop2(2));
objRaw_cropped = objRaw(REFcrop2(3):REFcrop2(4), REFcrop2(1):REFcrop2(2));

RefMultiply = 1;
if RefMultiply == 1
    objRawRow = reshape(objRaw, 1, numel(objRaw));
    RefPSF = load('/home/alingold/matlab/PSF_ref_20230821.mat');
    refPSF = RefPSF.refPSF;

    % Initialize objM as a column vector of zeros
    objM = zeros(size(refPSF, 1), 1);
    
    % Vectorized addition and scaling
    objM = refPSF * objRawRow(:);  % objRaw(:) ensures it's a column vector
    
    % Reshape the result into imageSize x imageSize
    objM = reshape(objM, [imageSize, imageSize]);
    % Plotting


    objMcrop = [95 170 105 180];
    objM = objM(objMcrop(3):objMcrop(4), objMcrop(1):objMcrop(2));
    objM = rot90(objM);
    % figure()
    % imagesc(objM);
    % colormap gray;
    % title('Multiplied by Ref');

end




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
% imagesc(objRaw);
% colormap gray;
% title('reconstruction');
imagesc(objM);
colormap gray;
title('Multiplied by Ref');

subplot(2,2,4);
imagesc(objRaw_cropped);
colormap gray;
title('reconstruction cropped to FOV');

function patches = image_to_patches(image, patchSize)
    [rows, cols] = size(image);
    patches = cell(ceil(rows/patchSize), ceil(cols/patchSize));
    for i = 1:patchSize:rows
        for j = 1:patchSize:cols
            patches{ceil(i/patchSize), ceil(j/patchSize)} = image(i:min(i+patchSize-1,rows), j:min(j+patchSize-1,cols));
        end
    end
end

function reconstructedImage = patches_to_image(patches)
    reconstructedImage = cell2mat(patches);
end
