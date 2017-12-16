function [meanLum, rmsContrast, weightedAverageMask, varargout] = getStatForPatch(varargin)
% Computes the luminance, contrast of a given image
im = double(varargin{1});
if nargin == 2
    fname = varargin{2};
end

[sizeX, sizeY] = size(im);

% Mean Lum
meanLum = nanmean(im(:));

% Contrast
normIm = (im./255.);
rmsContrast = nanstd(normIm(:));


% Centroid Spatial Frequency
intIm = im;
intIm(logical(isnan(im))) = 128;
Fim = fftshift(fft2(double(intIm)-meanLum));
nx = size(Fim, 2);
ny = size(Fim, 1);
Fmag = (abs(Fim)).^2;
f=-nx/2:nx/2-1;
[X Y]=meshgrid(f,f); % get coordinates of the power spectrum image
[theta rho]=cart2pol(X,Y); % equivalent in polar coordinates
rho=round(rho);
f1=zeros(nx/2+1,1);
fF=0:nx/2;
i = {};
for r=0:nx/2
    i{r+1}=find(rho==r); % for each freq return the location in the polar
    f1(r+1)=mean(Fmag(i{r+1})); % average power values of all the same freq
end
% convert the frequencies to cycles per degree
fCyclePerDegree = fF/(distDegrees(1)*size(im,2));
sel = fCyclePerDegree>0;
weights = f1(sel);
weightedAverageMask = sum(fCyclePerDegree(sel)' .* weights) ./ (sum(weights)+1e-10);

% Homogeneity
if nargout == 4
    patchHet = getHeterogeneityCagli(im, fname);
    varargout{1} = patchHet;
end

end