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

load(['C:\Cannula Microscope\IVCCM Software V3\Matlab Library\20230901_PSF.mat']);
RefPSF = load('C:\Cannula Microscope\IVCCM Software V3\Matlab Library\20230901_PSF_ref.mat');
%load(['20230731.mat']);
% 20230812
%%
tic;
%close all
%load test image CCM raw
%load(['/media/alingold/MenonLab/20230807_50um_25cal_100umWide/x24_y24_z0.mat'],'objCCM', 'objRef');
%load(["/media/alingold/MenonLab/20230822_5x5_multiBeadTest/x4_y0_z0.mat"],'objCCM', 'objRef');
%load(["/media/alingold/MenonLab/20230901_plant_Exposure1-3000ms/250ms.mat"],'objCCM', 'objRef');
%load(["C:\Cannula Microscope\20230822_plantTest\9.mat"],'objCCM', 'objRef');
%load(["/media/alingold/MenonLab/20230829_cottonFIber/x0_y100_z0.mat"],'objCCM', 'objRef');
%load(['/media/alingold/MenonLab/20230803_50um_bead_calTest/x1_y0_z0.mat'],'objCCM', 'objRef');
%load(['/media/alingold/MenonLab/20230731_50um_25cal/x20_y10_z0.mat'],'objCCM', 'objRef');
%load(['/media/alingold/MenonLab/CCM Ruipeng/2022_02_21_psf_for_plant_1/x15_y15_z0.mat'],'imgHDR', 'objRef');
load(["C:\Cannula Microscope\20230822_5x5_multiBeadTest\x4_y0_z0.mat"],'objCCM', 'objRef');
%load(["C:\Cannula Microscope\20230901_plant_Exposure1-3000ms/3000ms.mat"],'objCCM', 'objRef'); %1000ms for DC
%load(["C:\Cannula Microscope\20230825_plant2/4409.mat"],'objCCM', 'objRef'); 
%print(size(objCCM))

lambda = 4; %5 for no DC % 100000 for DC plant %10 for miltibead %50 for cotton
k = 200;
thresholdR = 0.001; % DC0.0005; % for thresholdR for 20230901_PSF 0.0012;
gamma_value = .8;
% 20230713 REFcrop = [400, 840, 440, 860];
% 20230821 REFcrop = [420, 1050, 220, 900]; 
REFcrop = [ 240, 960, 200, 940]; %crop reference to size of PSF
% 821 REFcrop2 = [35 90 45 100]; 
REFcrop2 = [ 40, 85, 23, 68]; % crop reference to the size of cannula

% fix broken pixel
objCCM(313, 457) = mean(mean(objCCM));

% find max and min within mask
masked_values = objCCM .* mask;
min_value_within_mask = min(min(masked_values));
max_value_within_mask = max(max(masked_values));

objCCM = (objCCM - min_value_within_mask) / (max_value_within_mask - min_value_within_mask);% Normalizing the image
objCCM = objCCM .* mask; % Applying the mask
objCCM = objCCM(top_left_row:bottom_right_row, top_left_col:bottom_right_col); % crop to size
objCCM = imresize(objCCM, [imageSize imageSize]); % resize to imsize
objCCM = objCCM - DC; % Subtract DC

objCCM = single(objCCM) - DC; % Convert objCCM to single-precision and subtract DC (likely a constant offset or bias).
frames.imgPro = reshape(objCCM, [1, imageSize^2])'; % Reshape objCCM into a column vector and transpose it. Store in frames.imgPro.
minImg = min(frames.imgPro(imgFOV)); % Find the minimum value in frames.imgPro for the region specified by imgFOV (likely a logical mask)
frames.imgPro(~imgFOV) = minImg; % Set the values in frames.imgPro where imgFOV is false to minImg, effectively masking out those regions.
y = double(frames.imgPro); % Convert frames.imgPro to double-precision for further mathematical operations and store in variable y.

clear refPro tform minImg; % Clear variables refPro, tform, and minImg from the workspace to free up memory.

if (min(lambda)<0) % Check if the minimum value in array lambda is negative.
  error('Illegal regularization parameter lambda') % If so, throw an error indicating an illegal regularization parameter.
end

% tikhonov
n = size(V,1); % Get the number of rows in matrix V and store it in variable n.
beta = U'*y; % Multiply the transpose of matrix U with vector y and store the result in variable beta.
zeta = s(:,1).*beta; % Element-wise multiply the first column of matrix s with vector beta and store in variable zeta.
ll = length(lambda); % Get the length of the lambda array and store it in variable ll.
x_lambda = zeros(n,ll); % Initialize a zero matrix with dimensions [n, ll] and store it in variable x_lambda.

    
% Treat each lambda separately.
for i=1:ll
    x_lambda(:,i) = V*(zeta./(s.^2 + lambda(i)^2)); 
    %x_lambda(:,i) = V*(zeta./(s + lambda(i))); L1 regularization is much worse
end

% extract guess from x_lambda
objRaw = x_lambda;
objRaw = reshape(objRaw,[psfSize psfSize zSize]);
% objRaw = abs(objRaw);
toc

L0 = 0; %L0 sparsity
if L0 == 1;
    tic
    % Initialize the parameters
    Omega = []; % Index set to keep track of selected columns of A
    r = y;      % Initial residual
         % Maximum number of non-zero coefficients or stopping criterion
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
end

% set threshold for objRaw
for n = 1:size(objRaw)
    for m = 1:size(objRaw)
        if objRaw(n, m) < thresholdR
            objRaw(n, m) = 0;
        end
    end
end


Ref = objRef; % pull ref img (ground truth)
objRef = objRef(REFcrop(3):REFcrop(4), REFcrop(1):REFcrop(2)); % crop to PSF
objRef = rot90(objRef); % rotate
objRef = imresize(objRef, [psfSize, psfSize]); % resize to psfSize
objRef = objRef(REFcrop2(3):REFcrop2(4), REFcrop2(1):REFcrop2(2)); % crop to illumminated region
objRaw_cropped = objRaw(REFcrop2(3):REFcrop2(4), REFcrop2(1):REFcrop2(2)); % crop reconstruction

RefMultiply = 1;
if RefMultiply == 1
    objRawRow = reshape(objRaw, 1, numel(objRaw)); % reshape objRaw into vector
    
    refPSF = RefPSF.refPSF; % extract refPSF

    objM = zeros(size(refPSF, 1), 1); % Initialize objM as a column vector of zeros
    objM = refPSF * objRawRow(:);  % objRaw(:) ensures it's a column vector
    objM = reshape(objM, [imageSize, imageSize]); % Reshape the result into imageSize x imageSize

    objMcrop = [95 170 105 180]; % crop reconstruction multiplied by ref PSF
    objM = objM(objMcrop(3):objMcrop(4), objMcrop(1):objMcrop(2));
    objM = rot90(objM);
    objM = abs(objM);

% apply gamma correction
for n = 1:size(objM, 1)
    for m = 1:size(objM, 2)
        objM(n, m) = objM(n, m) ^ (1 / gamma_value);
    end
end

else 
    objM = zeros(size(objRaw_cropped, 1), 1); % of RefMultipy == 1 set objM to zeros
end

toc

figure;
subplot(2,2,1);
imagesc(objRef);
colormap gray;
title('cropped ground truth to FOV');

subplot(2,2,2);
imagesc(objCCM);
colormap gray;
title('Processed Fiber');

subplot(2,2,3);
imagesc(objM);
colormap gray;
title('Multiplied by Ref');

subplot(2,2,4);
imagesc(objRaw_cropped);
colormap gray;
title(['reconstruction cropped to FOV lambda: ' num2str(lambda)]);
impixelinfo