%% Variable Needed
% Calibration data
% sfX : Image shift along X
% sfY : Image shift along Y
% imgRaw : HDR cannula image

frames(f).imgPro = frames(f).imgRaw;
if (frames(f).sfX ~= 0 || frames(f).sfY ~= 0)
    refPro = imref2d(size(frames(f).imgPro));
    tform  = affine2d([1 0 0; 0 1 0; frames(f).sfX frames(f).sfY 1]);
    frames(f).imgPro = imwarp(frames(f).imgPro,tform,'OutputView',refPro,'FillValues',0);
end
frames(f).imgPro = imresize(frames(f).imgPro,size(calibration.imgFOV));

if size(frames(f).imgPro) ~= size(calibration.imgFOV)
    error('Image size does not match with calibration. Check your image subregion setting.');
end

frames(f).imgPro = frames(f).imgPro - calibration.DC;
minImg = min(frames(f).imgPro(calibration.imgFOV));
frames(f).imgPro(~calibration.imgFOV) = minImg;
frames(f).imgPro = double(frames(f).imgPro);

clear refPro tform minImg;