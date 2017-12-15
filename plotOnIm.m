function plotOnIm(allCfg, allFiles, thisFile)
RFs = allFiles(1).RFs;
% im = imread(allFiles(1).imName);
% plot stuff on the image
if ~isfield(allCfg, 'print'); allCfg.print = false; end;
if ~isfield(allCfg, 'outputfile'); allCfg.outputfile = '.'; end;
savefile = allCfg.outputfile;
% make a screen for the borders
fullScreen = ones(1050, 1680)*128;
screenSize = [1680 1050];
fixPoint = screenSize/2;
% stimCenter = screenSize/2+52*[2 3];
stimCenter = allFiles(1).stimCenter;
% imSize = size(im);
nChan = size(RFs, 2);
nCond = length(allFiles);

% now the plots
upRatio = 10;
plotSize = [10 10]*upRatio;

if strcmp(allCfg.layout, 'channels')
    thisBaseline = allFiles(1).powerBaseline.powspctrm;
    % find the gray image
    allNormalized = [];
    data = cat(2, allFiles.lfpPower);
    data = cat(3, data.powspctrm);
    for cnd =1:nCond
        allNormalized = [allNormalized mean(data(:, :, cnd)./thisBaseline, 2)];
    end
    
    [~, b] = min(allNormalized, [], 2); probGrayInd = mode(b); %% Dangerous but logical
    thisBaseline = data(:,:,probGrayInd);
    
    for pl=1:nCond
        h = figure; if allCfg.print; set(h, 'visible', 'off'); end;
        imFull = imread(allFiles(pl).imName);
        fullScreen = imresize(imFull, size(imFull)*upRatio);
        imSize = size(fullScreen);
        imagesc(fullScreen); hold on; colormap gray; axis image;
        if strcmp(thisFile, 'lfpPower')
            thisPlot = allFiles(pl).lfpPower.powspctrm;
            data = log10(thisPlot./(thisBaseline));
            tfreq = allFiles(1).lfpPower.freq;
            tlabel = allFiles(1).lfpPower.label;
            for ch=1:nChan
                if ~(strcmp(tlabel(ch), 'V1-X') || strcmp(tlabel(ch), 'V1-84'))
                    plot(upRatio*(RFs(ch).centerposx-stimCenter(1))+imSize(1)/2 , ...
                        upRatio*(RFs(ch).centerposy-stimCenter(2))+imSize(1)/2 , 'g.')
                    xBorder = upRatio*(RFs(ch).centerposx - stimCenter(1)) +imSize(1)/2 - plotSize(1)/2 + 1;
                    yBorder = upRatio*(RFs(ch).centerposy - stimCenter(2)) +imSize(2)/2 + plotSize(2)/2 ;
                    tfr = tfreq/max(tfreq) * plotSize(1) + xBorder;
                    sel = 10 < tfreq & tfreq < 100;
                    tData = yBorder - data(ch, sel)/max(max(data(:, sel))+1e-20) * plotSize(2);
                    plot(tfr(sel), tData, 'r','LineWidth', 1); hold on;
                end
            end
        elseif strcmp(thisFile, 'stSpec')
            thisBaseline = allFiles(pl).powerBaseline.powspctrm;
            thisPlot = allFiles(pl).lfpPower.powspctrm;
            if isfield(allFiles, 'lfpPower2')
%                 powerLFP2 = cat(2, allFiles.lfpPower2);
%                 data2 = cat(3, powerLFP2.powspctrm);
                data2 = allFiles(pl).lfpPower2.powspctrm;
                thisPlot = (thisPlot+ data2)/2;
                disp('averaging two windows!!!');
            end
            data = log10(thisPlot./(1e-20+thisBaseline));
            tfreq = allFiles(1).lfpPower.freq;
            tlabel = allFiles(1).lfpPower.label;
            for ch=1:nChan
                if ~(strcmp(tlabel(ch), 'V1-X') || strcmp(tlabel(ch), 'V1-84'))
                    thisChan = str2num(tlabel{ch});
%                     nd = find((layout==thisChan));
%                     neighbors = layout(neighbourND(find(layout==thisChan), size(layout), [1 1]));
%                     nind = find(ismember(cellfun(@(x) str2num(x), tlabel), neighbors));
                    data = allFiles(ch, cnd);
                    data = cat(3, data.ppc1);
                    %                     plot(RFs(ch).centerposx, RFs(ch).centerposy, 'g.')
                    xBorder = RFs(ch).centerposx - stimCenter(1) +imSize(1)/2 - plotSize(1)/2 + 1;
                    yBorder = RFs(ch).centerposy - stimCenter(2) +imSize(2)/2 + plotSize(2)/2 ;
                    tfr = tfreq/max(tfreq) * plotSize(1) + xBorder;
                    tData = yBorder - data(ch, :)/(max(data(ch, :))+1e-20) * plotSize(2);
                    plot(tfr, tData, 'r','LineWidth', 2); hold on;
                end
            end
        end
        % zoom in to the plots
        borders = [[min(upRatio*([RFs.centerposx] - stimCenter(1)) +imSize(1)/2) ...
            min(upRatio*([RFs.centerposy] - stimCenter(2)) +imSize(2)/2)]-30*upRatio; ...
            [max(upRatio*([RFs.centerposx] - stimCenter(1)) +imSize(1)/2) ...
            max(upRatio*([RFs.centerposy] - stimCenter(2)) +imSize(2)/2)]+30*upRatio];
        xlim([borders(1,1) borders(2, 1)]);
        ylim([borders(1,2) borders(2, 2)]);
        
        % save fig
        figname = sprintf('cond%02d_lfpPowerOnIm.png', pl);
        print(h, fullfile(savefile, figname), '-dpng', '-r300');
    end
else
    warning('check layout parameter')
end