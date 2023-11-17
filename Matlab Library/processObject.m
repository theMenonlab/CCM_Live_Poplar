%% Variable Needed
% Calibration data
% gain
% threshold
% brightness
% sigma

% Thresholding and etc.
objProScale = frames(f).objRaw;
objProScale = objProScale ./ (max(objProScale(:)));
objProScale = objProScale*frames(f).gain;
objProScale(objProScale<frames(f).threshold) = frames(f).threshold;
objProScale = objProScale + frames(f).brightness;
objProScale(objProScale>1) = 1;

zIT =  1:calibration.zSize;
clear objProTempZ;
for z = 1:calibration.zSize
    objProTemp = objProScale(:,:,z);
    objSize2D = [size(frames(f).objRaw,1) size(frames(f).objRaw,2)];
    
    % Ext Compenstation
    %if frame.exCom == 1
    %    objProTemp = objProTemp./flipud(normArray);
    %    objProTemp = objProTemp./max(objProTemp(:));
    %    objProTemp = double(objProTemp);
    %end

    % Boundary masking
    if (frames(f).maskOpt == 1)
        padSize  = 100-frames(f).maskSize;
        objMask = padarray(fspecial('disk',frames(f).maskSize'),[padSize padSize]);
        
        objMask = imresize(objMask,objSize2D);
        objMask = circshift(objMask,[frames(f).maskShiftY,frames(f).maskShiftX]);
        objMask = (objMask == 0);
        objMask = ~objMask;
        objProTemp(~objMask) = 0;
    end

    % Convolution --  Need Fix
    if (frames(f).convOpt == 1)
        stepSize = 2;                 % Calibration step size 2um
        m = 5;                          % Oversampling factor
        pixelSize = stepSize/m;         % Oversampled step size
        sigma = 1*stepSize/pixelSize;   % 4um diameter bead
        
        % Disk convolution
        objOver = imresize(objProTemp,m*objSize2D);
        h = fspecial('disk', sigma);
        objProTemp = conv2(objOver,h,'same');
        objProTemp(objProTemp<frames(f).threshold) = frames(f).threshold;
    end
    objProTempZ(:,:,z) = flipud(objProTemp);
end
frames(f).objPro = objProTempZ;

clear objProTemp objProTempZ objMask objProScale padSize h objOver  ... 
      sigma pixelSize stepSize m