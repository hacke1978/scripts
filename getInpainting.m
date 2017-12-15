function out = getInpainting(allCfg)
% run inpainting for im
cx = allCfg.centerx;
cy = allCfg.centery;
maskSize = allCfg.cSize;
surroundSize = allCfg.sSize;
ssimSize = allCfg.ssimSize;
testImageName = allCfg.inputfile;

% fixed params
psz = 9;

% out the params
out.cx = cx;
out.cy = cy;
out.maskSize = maskSize;
out.surroundSize = surroundSize;
out.ssimSize = ssimSize;
out.inputfile = allCfg.inputfile;
out.name = allCfg.name;
out.originalName = allCfg.originalName;
out.stimName = allCfg.stimName;
out.name = allCfg.name;
out.ch = allCfg.ch;
out.statType = allCfg.statType;
out.psz = psz;

% run the algorithm
origImg = imread(testImageName);
imageSize = [600 600];
% crop original image to stimuli size
[sizeVal, sizeInd] = min(size(origImg(:, :, 1)));

% crop the image
if sizeVal > 1600
    origImg = imresize(origImg, 0.5);
    origImg = origImg(round(size(origImg, 1)/2-imageSize(1)/2:size(origImg, 1)/2+imageSize(1)/2-1), ...
        round(size(origImg, 2)/2-100-imageSize(2)/2:size(origImg, 2)/2-100+imageSize(2)/2-1), :);
elseif 802 < sizeVal && sizeVal < 1400
    origImg = origImg(round(size(origImg, 1)/2-imageSize(1)/2:size(origImg, 1)/2+imageSize(1)/2-1), ...
        round(size(origImg, 2)/2-100-imageSize(2)/2:size(origImg, 2)/2-100+imageSize(2)/2-1), :);
else
    origImg = origImg(round(size(origImg, 1)/2-imageSize(1)/2:size(origImg, 1)/2+imageSize(1)/2-1), ...
        round(size(origImg, 2)/2-imageSize(2)/2:size(origImg, 2)/2+imageSize(2)/2-1), :);
end

imCrop = zeros(maskSize, maskSize, 3);
imCrop = origImg(cy-surroundSize/2+1:cy+surroundSize/2, ...
    cx-surroundSize/2+1:cx+surroundSize/2, :);
mask = zeros(surroundSize, surroundSize);
mask(surroundSize/2-maskSize/2+1:surroundSize/2+maskSize/2, ...
    surroundSize/2-maskSize/2+1:surroundSize/2+maskSize/2) = 1;

if allCfg.make_movie
    [inpaintedImg,c,d,fillingMovie] = inpainting(imCrop,mask,psz);
    out.movie = fillingMovie;
else
    [inpaintedImg, ~, ~, ~] = inpainting(imCrop,mask,psz);
end

% ssim works on b&w images
grayInpainted = rgb2gray(double(inpaintedImg)/255);
grayOrig = rgb2gray(double(imCrop)/255);

% crop to ssim comparison size
cropInp = grayInpainted(surroundSize/2-ssimSize/2+1:surroundSize/2+ssimSize/2, ...
    surroundSize/2-ssimSize/2+1:surroundSize/2+ssimSize/2);
cropGray = grayOrig(surroundSize/2-ssimSize/2+1:surroundSize/2+ssimSize/2, ...
    surroundSize/2-ssimSize/2+1:surroundSize/2+ssimSize/2);

[fsim, fmask] = ssim(cropInp, cropGray);
out.fsim = fsim;
out.fmask = fmask;

if allCfg.make_image
    out.image_inpaint = inpaintedImg;
    out.image_original = imCrop;
end
end

