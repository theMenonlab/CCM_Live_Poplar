
%% List of necessary calibration variables : 
% All variables       : A, DC, U, s, V, calID, imageSize, imgFOV, normArray, objFOV, path, psfSize, scaleFactor
% Most necessary oens : U, s, V, DC, scaleFactor, imgFOV, objFOV, normArray

% 2D cal files does not have z parameter
if exist('zSize','var')==0
    zSize = 1;
end

if exist('calibration','var')==0
    error('Calibration is not loaded properly');
end