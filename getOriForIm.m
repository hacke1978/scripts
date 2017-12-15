function [peakOri] = getOriForIm(im, rfMask)
% get the main orientation using gabor filters

%%Parameters
orientation = [0:11.25:179];
% orientation = [0:45:179];
% orientation = [0:22.5:179];

freq = [.5:.5:2]/degDistances(1); % cycles per pixels
wavelength = round(1./freq); % pixels per cycle

% g = gabor(wavelength, orientation, 'SpatialFrequencyBandwidth', [.1 8]);
g = gabor(wavelength, orientation);
I = im;

%Functions
[mag,phase] = imgaborfilt(I,g);
nanMag = NaN(size(mag));

repMask = repmat(rfMask, 1, 1, size(mag, 3));
nanMag(repMask) = mag(repMask);

pwr = zeros(length(orientation), 1); ip = 0;
for ii=1:length(freq):size(nanMag, 3)
    ip=ip+1;
    thisOri = nanMag(:,:,ii:ii+ length(freq)-1);
    pwr(ip) = nansum(thisOri(:));
end

[mval, mind] = max(pwr);
peakOri = orientation(mind);