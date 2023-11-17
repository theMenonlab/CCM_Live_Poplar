
%imwrite(frames.imgRaw, 'm0_15exp.jpg')

 x{1}=('mix0.15exp.mat');
 x{2}=('mix0.1exp.mat');
 x{3}=('mix0.05exp.mat');
 
 ex=[0.15, 0.1, 0.05];
 
for i=1:3
    load(x{i})
    imgs=frames.imgRaw;
    
    rgbImg=zeros(size(imgs,1),size(imgs,2),3);
    rgbImg(:,:,1)=imgs./max(imgs);
    rgbImg(:,:,2)=imgs./max(imgs);
    rgbImg(:,:,3)=imgs./max(imgs);
     

    imwrite(rgbImg, ['m',num2str(i),'exp.tiff']);
%      J = imadjust(imgs,[],[],0.9);
%      figure;
%      imshow(J)

    hdr = lognormal(imgs);
    Lab = sRGB2Lab(hdr);
    Luminance = Lab(:,:,1);
    Luminance = adapthisteq(Luminance, 'NumTiles', numtiles);
    Luminance = imadjust(Luminance, LRemap, [0 1]);
    Lab(:,:,2) = Lab(:,:,2) * saturation;
    Lab(:,:,3) = Lab(:,:,3) * saturation;
end
 
files={['m1exp.tiff'],['m2exp.tiff'],['m3exp.tiff']};
hdr=makehdr(files, 'RelativeExposure', ex ./ ex(1));
final=hdr(:,:,1);
save('hdr1.mat', 'final');
% ends=rgb2gray(imread('hdr')); 

rgb = tonemap(hdr);
figure
imshow(rgb)

