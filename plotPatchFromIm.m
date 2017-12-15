function plotPatchFromIm(allCfg, allFiles)
%Get imStats for a given image

screenSize = [1680 1050];
fixPoint   = screenSize/2;
RFs = allFiles(1).RFs;
imSize = size(imread(allFiles(1).imName));
imCenter = imSize/2;
nCond = length(allFiles);
nChan = length(RFs);
savename = allCfg.outputfile;
caccept = allFiles(1).caccept;

if isfield(allFiles(1), 'sortInd'), sortInd = allFiles(1).sortInd; else sortInd = 1:nCond; end;

if strcmp(allCfg.name, 'Hermes')
    layout = reshape(64:127, 8, 8);
    if strcmp(allCfg.type, 'grating-ori')
        stimCenter = round(fixPoint+51*([2 5]));
    else
%         stimCenter = round(fixPoint+degDistances(1)*([2 3]));
                  stimCenter = round(fixPoint+degDistances(1)*([2 6.5]));
    end
    label = [RFs.label];
    label = cellfun(@(x) str2num(x(4:end)), label);
elseif strcmp(allCfg.name, 'Isis')
    layout = vertcat([NaN, NaN, 27, 1, NaN, NaN],[reshape(17:26, 5, 2) [28:32]' reshape(2:16, 5, 3)]);
    stimCenter = [960+16 660-136];
    label = vertcat({RFs.label});
    label = cellfun(@(x) str2num(x(4:end)), label);
elseif strcmp(allCfg.name, 'Ares')
    layout = vertcat([NaN, NaN, 27, 1, NaN, NaN],[reshape(17:26, 5, 2) [28:32]' reshape(2:16, 5, 3)]);
    stimCenter = [980 520];
    label = vertcat({RFs.label});
    label = cellfun(@(x) str2num(x(4:end)), label);
end

% get max RF/ surround size
cSize = 2*round(nanmedian([RFs.sigmaX RFs.sigmaY]));
sSize = 2*cSize;

if strcmp(allCfg.layout, 'channels')
    nr = size(layout, 1); nc = size(layout, 2);
    for cnd = 1:nCond
        cnd = sortInd(cnd);
        fullScreen = 128*ones(screenSize(2), screenSize(1));
        im = imread(allFiles(cnd).imName);
        imSize = size(im);
        fullScreen(stimCenter(2)-imSize(2)/2+1:stimCenter(2)+imSize(2)/2, ...
            stimCenter(1)-imSize(1)/2+1:stimCenter(1)+imSize(1)/2, :) = im;
        patchPerCond = 128*ones(nr*sSize, nc*sSize, size(im, 3));
        for ch = 1:nChan
            [rind, cind] = ind2sub(size(layout), find((layout==label(ch))));
            % RFs w.r.t to fullScreen center
            thisPatch = getRFPatch(fullScreen, RFs, ch, sSize);
            patchPerCond((rind-1)*sSize+[1:sSize], (cind-1)*sSize+[1:sSize], :) = thisPatch;
        end
        h = figure(); set(h, 'visible', 'off');
        imagesc(patchPerCond); colormap gray; axis image; hold on;
        for ch = 1:nChan
            [rind, cind] = ind2sub(size(layout), find((layout==label(ch))'));
            ellipsedraw(RFs(ch).sigmaX, RFs(ch).sigmaY, (rind-1)*sSize+sSize/2, (cind-1)*sSize+sSize/2,-RFs(ch).angle, 'r', [sSize sSize], 0);
            hold on;
            text((rind-1)*sSize+sSize/2, (cind-1)*sSize+sSize/2, {label(ch)}, 'Color', 'r');
            set(gca, 'FontSize', 5);
        end
        figname = sprintf('cond%02d_imagePatch.png', cnd);
        print(h, fullfile(savename, figname), '-dpng', '-r300');
    end
elseif strcmp(allCfg.layout, 'stimuli')
    [im, ~, alphaMask] = imread(allFiles(1).imName);
    im(repmat(alphaMask, [1 1 3])==0) = 128;
    if size(im, 3) == 3
        nr = ceil(sqrt(nCond)); nc = ceil(sqrt(nCond));
    else
        nr = ceil(sqrt(nCond)); nc = ceil(sqrt(nCond));
    end
    patchPerCha = 128*ones(nr*sSize, nc*sSize, size(im, 3), nChan);
    rind = 1; cind = 0;
    for cnd = 1:nCond
        cnd = sortInd(cnd);
        if cind == nc
            cind = 0;
            rind = rind + 1;
        end
        
        cind = cind + 1;
        fullScreen = 128*ones(screenSize(2), screenSize(1), size(im, 3));
        [im, ~, alphaMask] = imread(allFiles(cnd).imName);
        imSize = size(im);
        im(repmat(alphaMask, [1 1 3]) == 0) = 128;
        if size(im, 3) == 1
            im = repmat(im, 1, 1, 3);
        end
        fullScreen(stimCenter(2)-imSize(2)/2+1:stimCenter(2)+imSize(2)/2, ...
            stimCenter(1)-imSize(1)/2+1:stimCenter(1)+imSize(1)/2, :) = im;
        for ch = 1:nChan
            if caccept(ch)
                % RFs w.r.t to fullScreen center
                thisPatch = getRFPatch(fullScreen, RFs, ch, sSize);
                patchPerCha((rind-1)*sSize+[1:sSize], (cind-1)*sSize+[1:sSize], :, ch) = thisPatch;
            end
        end
    end
    for ch = 1:nChan
        if caccept(ch)
            h = figure(); set(h, 'visible', 'off');
            if size(im, 3) == 3
                imagesc(patchPerCha(:, :, :, ch)./255);
            else
                imagesc(patchPerCha(:, :, ch))
                colormap gray;
            end
            axis image; hold on;
            rind = 1; cind = 0;
            for cnd = 1:nCond
                cnd = sortInd(cnd);
                if cind == nc
                    cind = 0;
                    rind = rind + 1;
                end
                cind = cind + 1;
                ellipsedraw(RFs(ch).sigmaY, RFs(ch).sigmaX, (cind-1)*sSize+sSize/2, (rind-1)*sSize+sSize/2,-RFs(ch).angle, 'r', [sSize sSize], 0);
                hold on;
                text((cind-1)*sSize+sSize/2, (rind-1)*sSize+sSize/2, {label(ch)}, 'Color', 'r');
                set(gca, 'FontSize', 5);
            end
            figname = sprintf('ch%02d_imagePatch.png', label(ch));
            print(h, fullfile(savename, figname), '-dpng', '-r300');
        end
    end
end