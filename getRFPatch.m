function [rfPatch, pMask, XX, YY] = getRFPatch(im, RFs, ch, sSize)
enlargenScreen = [2000 2000]; % temporary fix for insane RFs
screenSize = [1680 1050]+enlargenScreen;
fullScreen = 128*ones(screenSize(2), screenSize(1), size(im, 3));
imSize = size(squeeze(im(:, :, 1)));
imCenter = screenSize/2;
fullScreen(imCenter(2)-imSize(1)/2+1:imCenter(2)+imSize(1)/2, ...
            imCenter(1)-imSize(2)/2+1:imCenter(1)+imSize(2)/2, :) = im;
fixPoint   = screenSize/2;
% stimCenter = [960+16 660-136]; 

% ty = screenSize(2) - [RFs.centerposy];
% get RF positions wrt fixation point
fx = round([RFs.centerposx] + enlargenScreen(1)/2) - fixPoint(1);
fy = round([RFs.centerposy] + enlargenScreen(2)/2) - fixPoint(2);

% % get RF positions wrt stimcenter
% stimShift = stimCenter - fixPoint;
% sx = round(fx - stimShift(1));
% sy = round(fy - stimShift(2));

% get max RF/ surround size
cSize = 2*round(median([RFs(ch).sigmaX RFs(ch).sigmaY]));
% sSize = 80;

rfPatch = fullScreen(imCenter(2) + fy(ch)-sSize/2+1:imCenter(2) + fy(ch)+sSize/2, ...
                   imCenter(1) + fx(ch)-sSize/2+1:imCenter(1) + fx(ch)+sSize/2, :);
% figure, imagesc(rfPatch); set
% [pMask, XX, YY] = ellipsedrawMore(RFs(ch).sigmaX, RFs(ch).sigmaY, sSize/2, sSize/2,-RFs(ch).angle, 'r', [sSize sSize], 0);