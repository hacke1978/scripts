function getStatsPerChannel()
im = imread(imList{1});
[rfPatch, pMask, XX, YY] = getRFPatch(im, RFs, ch, stimCenter);
stats = getStatForIm(rfPatch);