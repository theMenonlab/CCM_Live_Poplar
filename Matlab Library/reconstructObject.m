%% Variable Needed
% Calibration data
% ctrlFactor : optimizatino parameter

y = double(frames(f).imgPro(calibration.imgFOV));
if frames(f).method == 0
    [x] = tikhonovFast(calibration.U,calibration.s,calibration.V,y,frames(f).ctrlFactor);
elseif frames(f).method == 1
    [x] = dbs(calibration.A, calibration.y, frames(f).ctrlFactor, 0);
elseif frames(f).method == 2
    [x] = tsvdFast(calibration.U,calibration.s,calibration.V,y,round(frames(f).ctrlFactor*length(calibration.s)/100));
elseif frames(f).method == 3
    [x imageNos2D] = dbsNoise(calibration.A, y, frames(f).ctrlFactor, 0, calibration.imgFOV, calibration.objFOVInd);
else
    [x] = tikhonovFast(calibration.U,calibration.s,calibration.V,y,frames(f).ctrlFactor);
end

frames(f).objRaw = zeros(1,calibration.psfSize^2*calibration.zSize);
frames(f).objRaw(calibration.objFOVInd) = x;
frames(f).objRaw = reshape(frames(f).objRaw,[calibration.psfSize calibration.psfSize calibration.zSize]);
frames(f).objRaw = abs(frames(f).objRaw);
