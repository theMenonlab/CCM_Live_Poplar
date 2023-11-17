 %% Init variables
clear; clc;
imageSize = 256;
zSize = 1;
sfX = [0]; % need to put in if there is a shift in callibration (oposite sign of labview)
sfY = [0]; % need to put in if there is a shift in callibration (oposite sign of labview)
psfSize = 100; %201 Ruipeng was doing more beads 200x200imgs 200um FOV 
scaleFactor = .5;

path{1} = 'C:\Cannula Microscope\20230901_PSF\';  zInd(1) = 0; % don't forget to change name at botom if needed
% load('Mask.mat'); % for normalization mask for ccm
% load('ObjMask.mat'); % mask for the beads 

filename0 = ['x50_y50_z0.mat'];
data0 = load([path{1},filename0]);

fprintf('Structure %d:\n', i);
fields = fieldnames(data0);
for j = 1:length(fields)
    fieldName = fields{j};
    fprintf('\t%s: %s\n', fieldName);
end
disp(['CCM size: ' num2str(size(data0.objCCM))]);
disp(['Ref size: ' num2str(size(data0.objCCM))]);

DC  = 0; % black point
% load('imgFOV');
%% create mask
% left_x = 122;
% right_x = 252;
% y_top = 118;
% 20230807
left_x = 232;
right_x = 470;
y_top = 225;

imgHDR = data0.objCCM;

% Resizing the image
%imgHDR = imresize(imgHDR, scaleFactor);

% Displaying the original image
figure;
%imshow(imgHDR/max(max(imgHDR)));
imagesc(imgHDR);

title('Original Image');
impixelinfo;

[rows, cols] = size(imgHDR);

% Calculating the center and diameter of the circular mask
center_x = (left_x + right_x) / 2;
mask_diameter = right_x - left_x;
center_y = y_top + mask_diameter/2;
radius = mask_diameter / 2;

% Creating the mask using logical indexing
[x, y] = meshgrid(1:cols, 1:rows);
mask = (x - center_x).^2 + (y - center_y).^2 <= radius^2;

%fix broken pixel
imgHDR(313, 457) = 0;

% Apply the mask and find the masked values
masked_values = imgHDR .* mask;

% Find the min and max value within the mask
min_value_within_mask = min(masked_values(masked_values ~= 0))
max_value_within_mask = max(max(masked_values))

% Find Second highest value
% A_reshaped = reshape(masked_values, [], 1);
% sorted_values = sort(A_reshaped, 'descend');
% max_value_within_mask = sorted_values(2);

% Normalizing the image
imgHDR = (imgHDR - min_value_within_mask) / (max_value_within_mask - min_value_within_mask);

% Applying the mask
masked_image = imgHDR .* mask;
%masked_image = masked_image - mean(mean(masked_image));
% Displaying the masked image
figure
imshow(imgHDR)
title('Masked Image')
impixelinfo
save('Mask.mat', 'mask')
%% crop to mask
% Finding the rows and columns that correspond to the mask
[row_mask, col_mask] = find(mask);

% Finding the top-left and bottom-right coordinates of the masked region
top_left_row = min(row_mask);
bottom_right_row = max(row_mask);
top_left_col = min(col_mask);
bottom_right_col = max(col_mask);

% Cropping the image to the bounding box of the mask
cropped_image = masked_image(top_left_row:bottom_right_row, top_left_col:bottom_right_col);
mask_cropped = mask(top_left_row:bottom_right_row, top_left_col:bottom_right_col);
% Displaying the cropped image
figure;
imshow(cropped_image);
title('Cropped Masked Image');
impixelinfo;

resized_image = imresize(cropped_image, [imageSize imageSize]);
mask_cropped_resized = imresize(mask_cropped, [imageSize imageSize]);
imgFOV = mask_cropped_resized;

figure
imshow(resized_image)
title('resized Masked Image')
impixelinfo

%% subtract background

load(["C:\Cannula Microscope\20230901_background\x0_y0_z0.mat"],'objCCM');

background = objCCM;

figure()
imshow(background / max(max(background)))
title('background')
impixelinfo
%% Get Calibration CCM
psf = zeros(imageSize,imageSize,psfSize,psfSize,zSize,'single');

for z = 0:zSize-1 % z is 1 for me
    fprintf('-------------------------z = %d------------------------\n',z);
    for y = 0:psfSize-1
        y1 = y;
        fprintf('y = %d\n',y);
        for x = 0:psfSize-1
            x1 = x;
            filename = ['x',num2str(x1), '_y', num2str(y1), '_z', num2str(zInd(z+1))];
            data = load([path{z+1},filename]);
            
            % Print the fields and values for each element in the structure array

    
            imgHDR = data.objCCM;
            
            if isempty(imgHDR) 
                fprintf('Warning: imgHDR is empty for x=%d, y=%d, z=%d\n', x, y, z);
                continue; % Skip to the next iteration
            end
            
            %subtrack background
            %imgHDR = imgHDR - background;

            %fix broken pixel
            imgHDR(313, 457) = mean(mean(imgHDR));
            % Apply the mask and find the masked values
%             masked_values = imgHDR .* mask;
%             min_value_within_mask = min(min(masked_values));
%             max_value_within_mask = max(max(masked_values));
%             % Normalizing the image
%             imgHDR = (imgHDR - min_value_within_mask) / (max_value_within_mask - min_value_within_mask);
            % Applying the mask
            masked_image = imgHDR .* mask;
            % masked_image = masked_image - mean(mean(masked_image));
            % crop to size
            cropped_image = masked_image(top_left_row:bottom_right_row, top_left_col:bottom_right_col);
            % resize to imageSize
            resized_image = imresize(cropped_image, [imageSize imageSize]);
            psf(:,:,y+1,x+1,z+1) = (resized_image - DC); % 5D PSF, img_x, img_y, X, Y, Z
        end
    end
end

%clear filename Exp imgHDR zInd
fprintf('1. PSF collation complete.\n');
size(psf)
%save('PSF.mat', 'psf', '-v7.3')


%% Get Calibration Ref
psf = zeros(imageSize,imageSize,psfSize,psfSize,zSize,'single');

for z = 0:zSize-1 % z is 1 for me
    fprintf('-------------------------z = %d------------------------\n',z);
    for y = 0:psfSize-1
        y1 = y;
        fprintf('y = %d\n',y);
        for x = 0:psfSize-1
            x1 = x;
            filename = ['x',num2str(x1), '_y', num2str(y1), '_z', num2str(zInd(z+1))];
            data = load([path{z+1},filename]);
            
            % Print the fields and values for each element in the structure array

    
            imgHDR = data.objRef;
            if isempty(imgHDR)
                fprintf('Warning: imgHDR is empty for x=%d, y=%d, z=%d\n', x, y, z);
                continue; % Skip to the next iteration
            end



            min_value = min(min(imgHDR));
            max_value = max(max(imgHDR));
            % Normalizing the image
            imgHDR = (imgHDR - min_value) / (max_value - min_value);

            % resize 
            resized_image = imresize(imgHDR, [imageSize imageSize]);
            psf(:,:,y+1,x+1,z+1) = (resized_image - DC); % 5D PSF, img_x, img_y, X, Y, Z
        end
    end
end

% Al's transform PSF
[imgX, imgY, X, Y, Z] = size(psf);
refPSF = zeros(imgX * imgY, X * Y * Z);

for x = 1:X
    for y = 1:Y
        for z = 1:Z
            positionIndex = (z-1) * X * Y + (y-1) * X + x;
            img = reshape(psf(:,:,x,y,z), [imgX * imgY, 1]);
            refPSF(:, positionIndex) = img;
        end
    end
end


%clear filename Exp imgHDR zInd
fprintf('1. PSF collation complete.\n');
size(refPSF)

save('20230901_PSF_ref.mat', 'refPSF', '-v7.3')

%% Mean Pixel Value for checking psf and plotting brightness distrubution from cannula
for z = 0:zSize-1
    for x = 0:psfSize-1
        for y = 0:psfSize-1
            temp = double(psf(:,:,y+1,x+1,z+1));

            normArray(y+1,x+1,z+1) = mean(temp(:)); %mean pixel value
        end
    end
    figure; imagesc(normArray(:,:,z+1)); colorbar; axis square; axis xy;
    title(sprintf('Mean pixel value, z=%d',z), 'FontName', 'Arial', 'FontSize', 14);
    impixelinfo
end

%% ** Find DC (black point, dark content) from images taken outside FOV and Save **
% Measure DC from experimentally collected PSF.
counter = 0;
DCA = zeros(imageSize);
for z = 0:zSize-1
    for x = 0:psfSize-1
        for y = 0:psfSize-1
            if (normArray(y+1,x+1) < 2000)  %0.1 if the mean pixel value is less than 0.13, it is set to background
                DCA = DCA + psf(:,:,y+1,x+1,z+1); %DCA is the sum of all images
                counter = counter+1;
            end
        end
    end
end
DC = DCA / counter; % turns sum into ave
DC(isnan(DC)) = 0;
figure; imagesc(DC); colorbar; axis square; axis xy;
impixelinfo
clear counter DCA;

%% Subtract DC 
% subtrack ave value from psf
for z = 0:zSize-1
    for x = 0:psfSize-1
        for y = 0:psfSize-1
            temp = psf(:,:,y+1,x+1,z+1);
            psf(:,:,y+1,x+1,z+1) = temp - DC;
        end
    end
end
fprintf('2. DC Subtraction complete.\n');

%% Find Ave image
temp = zeros(imageSize);
tempC = 0;
for z = 0:zSize-1
    for x = 0:psfSize-1
        for y = 0:psfSize-1
            temp = temp + psf(:,:,y+1,x+1,z+1);
            tempC = tempC + 1;
        end
    end
end
temp = temp ./ tempC;
figure;imagesc(temp);colorbar;

%% ** Define cannula image FOV **
%sample = temp;
%imgFOV = (sample > 0.025);%0.1
imgFOV = mask_cropped_resized;  % sometime you can directly use the mask to define the FOV
sample(~imgFOV) = 0;
figure;
subplot(1,2,1); imagesc(sample);
subplot(1,2,2); imagesc(imgFOV);

%% Apply imgFOV apply FOV crop
for z = 0:zSize-1
    for x = 0:psfSize-1
        for y = 0:psfSize-1
            temp = psf(:,:,y+1,x+1,z+1);
            temp(~imgFOV) = 0;
            psf(:,:,y+1,x+1,z+1) = temp;
        end
    end
end
fprintf('3. Image FOV Complete.\n');

%% Norm Pixel Value
clear normArray;
for z = 0:zSize-1
    for x = 0:psfSize-1
        for y = 0:psfSize-1
            temp = psf(:,:,y+1,x+1,z+1);
            
%             factor = 1;
%             r2 = (x - 51)*(x-51) + (y - 51)*(y-51);
%             if  r2 < 40 * 40
%                 factor =  20;
%             end
            
            normArray(y+1,x+1,z+1) =  norm(temp(:)); % norm each image, used to find CCM FOV, not used normally
        end
    end
    figure; imagesc(normArray(:,:,z+1)); colorbar; axis square; axis xy;
    title(sprintf('Norm pixel value, z=%d',z), 'FontName', 'Arial', 'FontSize', 14);
end

%% ** Define object FOV **
sample = normArray;
objFOV = (sample > 1.5);%20
% objFOV = ObjMask;
sample(~objFOV) = 0;

figure;
subplot(1,2,1); imagesc(objFOV(:,:,1)); axis xy;
subplot(1,2,2); imagesc(sample); axis xy;
objFOVInd = find(objFOV==1);
fprintf('4. Object FOV Complete.\n');

%% Correlation coefficient, Diagonal scan. Optional
% used to evaluate the quality of the PSF
for x = round(psfSize/2)
    for y = 1:psfSize
        temp = psf(:,:,y,y);
        psfAcrossD(:,y) = temp(:);
        
    end
end

corrMapD = corrcoef(psfAcrossD);
invdig = corrMapD(1:psfSize,psfSize-(0:(psfSize-1)));

figure;imagesc(corrMapD); colorbar;
title('Correlation coefficient map of Diagonal line');

figure; plot(corrMapD(:,30));
title('Correlation line scan at 30th diagonal psf');

figure; plot(corrMapD(:,45));
title('Correlation line scan at 45th diagonal psf');

figure; plot(diag(invdig));
title('Inverse Diagonal plot');

clear corrMapD psfAcrossD invdig
%% Al's transform PSF
[imgX, imgY, X, Y, Z] = size(psf);
psf2D = zeros(imgX * imgY, X * Y * Z);

for x = 1:X
    for y = 1:Y
        for z = 1:Z
            positionIndex = (z-1) * X * Y + (y-1) * X + x;
            img = reshape(psf(:,:,x,y,z), [imgX * imgY, 1]);
            A(:, positionIndex) = img;
        end
    end
end


%% Transform from 5D to 2D array. Normalize the 1D image here.
for z = 1:zSize
    for x = 1:psfSize
        %fprintf('x = %d\n', x);
        for y = 1:psfSize
            temp = psf(:,:,y,x,z);
            temp = temp(imgFOV);
%            temp = temp ./ norm(temp);
            psf4D(:,y,x,z) = temp;
        end
    end
end

A = zeros(size(psf4D,1),sum(objFOV(:)),'single');
for x = 1:size(psf4D,1)
    %fprintf('x = %d\n', x);
    temp = squeeze(psf4D(x,:,:,:));
    temp = temp(objFOV);
    A(x,:) = temp(:);
end
clear x y psf4D temp;
fprintf('5. 4D PSF to 2D A Complete.\n');
% save('A.mat', 'A', '-v7.3')
% Remove 4D psf and Save for later use
clear psf;
%%
load('A.mat', 'A')
%% svd
% Basic analysis of A
tic;
    % same as [U,s,V] = csvd(A);
    [m,n] = size(A);
    if (m >= n)
      [U,s,V] = svd(full(A),0); s = diag(s);
    else
      [V,s,U] = svd(full(A)',0); s = diag(s);
    end
toc;

fprintf('Cond # : %d\n', s(1)/s(end));
figure; semilogy(s);

clear m n;
fprintf('6. SVD Complete.\n');

% Clear some stuff
clear corrMapD psfAcrossDx y m n temp invdig ans sample

% Save Calibration Data
calID = randi([0 100000], 1, 1);

save('20230901_PSF_DC.mat','-v7.3');
