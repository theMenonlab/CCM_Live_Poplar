%% Variable Needed
% Every variables to be saved
% saveMode   : single or series

if strcmp(saveMode, 'single')
    f = 1;
    clear frames;
elseif strcmp(saveMode, 'series')
    f = f + 1;
end

frames(f).exposure      = Exp;
frames(f).xPos          = xPos;
frames(f).yPos          = yPos;
frames(f).zPos          = zPos;
frames(f).timeStamp     = datestr(now);

frames(f).gain          = gain;
frames(f).threshold     = threshold;
frames(f).brightness    = brightness;
frames(f).convOpt       = convOpt;
frames(f).maskOpt       = maskOpt;
frames(f).maskSize      = maskSize;
frames(f).maskShiftX    = maskShiftX;
frames(f).maskShiftY    = maskShiftY;

frames(f).sfX           = sfX;
frames(f).sfY           = sfY;
frames(f).ctrlFactor    = ctrlFactor;
frames(f).method        = method;

frames(f).imgRaw        = imgRaw;
frames(f).objRef        = objRef;
frames(f).imgPro        = 0;
frames(f).objRaw        = 0;
frames(f).objPro        = 0;
