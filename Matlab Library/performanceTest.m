%%
%for i = 1:10
%    [x] = tsvdFast(U,s,V,y,10);
%    [x] = tikhonovFast(U,s,V,y,10);
%end

%% Overall Performance test

refRaw = 1;
imgRaw = imgRawFrames{1};
affMat = [1 0 0; 0 1 0; sfX sfY 1];
tform = affine2d(affMat);

for j = 1:10
    method = 2;
    ctrlFactor = 10;
    calibrationCheck;
    processAndorImage;
    reconstructImage;
    processReconImage;
    saveFrame;
end