function [out] = getOrientation(allCfg)
% get the main orientation using gabor filters
cx = allCfg.centerx;
cy = allCfg.centery;
maskSize = allCfg.cSize;
surroundSize = allCfg.sSize;
ssimSize = allCfg.ssimSize;
imName = allCfg.stimName;
stimDir = allCfg.stimDir;
imDir= fullfile(stimDir, sprintf('stimulus%03d.png', str2num(imName(10:end-4))));
im = imread(imDir);

% Parameters
orientation = [0:11.25:179];

if strcmp(allCfg.name, 'hermes')
    pixPerDeg = 51;
elseif strcmp(allCfg.name, 'isis')
    pixPerDeg = 40;
elseif strcmp(allCfg.name, 'ares')
    pixPerDeg = 40;
end
freq = [.5:.5:3]/pixPerDeg; % cycles per pixels
wavelength = round(1./freq); % pixels per cycle

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

%% Get the orientation
g = gabor(wavelength, orientation);

if ndims(im)==3
    I = rgb2gray(im);
    disp('image has to be gray')
else
    I = im;
end

imCrop = zeros(maskSize, maskSize, 3);
imCrop =  I(cy-maskSize/2+1:cy+maskSize/2, ...
    cx-maskSize/2+1:cx+maskSize/2, :);

% get gabor outputs
[mag, phase] = imgaborfilt(imCrop, g);
pwr = zeros(length(orientation), 1); ip = 0;
for ii=1:length(freq):size(mag, 3)
    ip=ip+1;
    thisOri = mag(:,:,ii:ii+ length(freq)-1);
    pwr(ip) = sum(thisOri(:));
end

[mval, mind] = max(pwr);
peakOri = orientation(mind);

out.peakOri = peakOri;