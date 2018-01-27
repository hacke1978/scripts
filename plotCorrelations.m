function plotCorrelations(allCfg, out)

if strcmp(allCfg.corrLayout, 'group')
    nameInd = cellfun(@(x) length(x), out, 'UniformOutput', false);
    nameInd = [nameInd{:}];
    
    for nind = unique(nameInd) %
        % dumb
        if nind == 22
            savefile = '/mnt/v7k/projects/MWNaturalPredict/Figures/NatImFix/Isis';
        elseif nind == 64
            savefile = '/mnt/v7k/projects/MWNaturalPredict/Figures/NatImFix/Hermes';
        end
        
        allOut = cat(1, out{nameInd == nind});
        label = [allOut(1, :).label];
        nChan = size(allOut, 2);
        nCond = size(allOut, 1);
        for ch = 1:nChan
            h = figure(); if allCfg.print; set(h, 'visible', 'off'); end;
            thisOut = cat(2, squeeze(allOut(:, ch)));
            thisCh = str2num(label{ch});
            allRates = horzcat(thisOut.rates);
            allPeaks = horzcat(thisOut.peaks);
            allStats = horzcat(thisOut.stats);
           if allCfg.isOri
            allOri = horzcat(thisOut.ori);
           end
            
            % RATES
            % get correlations for rates
            [sCorr, sVal] = corr(allStats', allRates', 'type', 'Spearman');
            allSpearRate = [sCorr, sVal];
            [pCorr, pVal] = corr(allStats', allRates', 'type', 'Pearson');
            allPearsonRate = [pCorr, pVal];
            nbin = 30;
            eqBins = quantile(allStats,nbin-1);
            eqBins = [0 eqBins 1];
            meanRates = []; meanStats = []; meanPeaks = [];
            for ii=1:length(eqBins)-1
                sel = eqBins(ii) <= allStats & allStats < eqBins(ii+1);
%                 sum(sel)
                meanRates(ii) = mean(allRates(sel));
                meanPeaks(ii) = mean(allPeaks(sel));
                meanStats(ii) = mean(allStats(sel));
            end
            plot(meanStats, meanRates, 'ko')
            title(sprintf('Pearson R %.2f, pVal %.d , Spearman R %.2f, pVal %.d ', ...
                allPearsonRate(1), allPearsonRate(2),...
                allSpearRate(1), allSpearRate(2)), 'fontsize', 6, 'fontweight', 'bold');
            set(gca, 'FontSize', 6);
            
            % save rate stuff
            if allCfg.do_lfpPower2
                figname = sprintf('ch%02d_corrRate2.png', thisCh);
            else
                figname = sprintf('ch%02d_corrRate.png', thisCh);
            end
            print(h, fullfile(savefile, figname), '-dpng', '-r300');
            close(h);
            
            % PEAKS
            % get correlations for peaks
            [sCorr, sVal] = corr(allStats', allPeaks', 'type', 'Spearman');
            allSpearRate = [sCorr, sVal];
            [pCorr, pVal] = corr(allStats', allPeaks', 'type', 'Pearson');
            allPearsonRate = [pCorr, pVal];
            
            % plot peaks
            h = figure(); if allCfg.print; set(h, 'visible', 'off'); end;
            plot(meanStats, meanPeaks, 'ko'); %ylim([-2 20]);
            title(sprintf('Pearson R %.2f, pVal %.d , Spearman R %.2f, pVal %.d ', ...
                allPearsonRate(1), allPearsonRate(2),...
                allSpearRate(1), allSpearRate(2)), 'fontsize', 6, 'fontweight', 'bold');
            set(gca, 'FontSize', 6);
            
            % save peak stuff
            if allCfg.do_lfpPower2
                figname = sprintf('ch%02d_corrPeak2.png', thisCh);
            else
                figname = sprintf('ch%02d_corrPeak.png', thisCh);
            end
            print(h, fullfile(savefile, figname), '-dpng', '-r300');
            close(h);
            
% %             X = [ [out{1}(8).stats]'];
% %             y = [out{1}(8).peaks];
% %             fitlm(X, y);
% %             
            % GLM stuff
%              X = [allStats; allOri]';
%               y = allPeaks';
%            [bPred, devPred, statsPred] = glmfit(allStats, y);
%             [bOri, devOri, statsOri]  = glmfit(allOri, y);
%             [bFull, devFull, statsFull]  = glmfit(X, y);
%             yhatPred = glmval(bPred, allStats, 'identity');
%             yhatOri = glmval(bOri, allOri, 'identity');
%             yhatFull = glmval(bFull, X, 'identity');
            % %             figure
% %             plot(X, y, '.'); hold on
% %             plot(X, yhat, 'r.')
            
        end
    end
elseif strcmp(allCfg.corrLayout, 'channels')
    savefile = allCfg(1).outputfile;
    nChan = length(out);
    label = [out.label];
    allRates = vertcat(out.rates);
    allPeaks = vertcat(out.peaks);
    allStats = vertcat(out.stats);
    for ch=1:nChan
        % get Ch
        thisCh = str2num(label{ch});
        % RATES
        h = figure(); if allCfg.print; set(h, 'visible', 'off'); end;
        plot(allStats(ch, :), allRates(ch, :), 'ko'); %ylim([-2 20]);
        title(sprintf('Pearson R %.2f, pVal %.d , Spearman R %.2f, pVal %.d ', ...
            out(ch).allPearsonRate(1) , out(ch).allPearsonRate(2),...
            out(ch).allSpearRate(1) , out(ch).allSpearRate(2)));
        
        % save rate stuff
        if allCfg.do_lfpPower2
            figname = sprintf('ch%02d_corrRate2.png', thisCh);
        else
            figname = sprintf('ch%02d_corrRate.png', thisCh);
        end
        print(h, fullfile(savefile, figname), '-dpng', '-r300');
        close(h);
        
        % PEAKS
        h = figure(); if allCfg.print; set(h, 'visible', 'off'); end;
        plot(allStats(ch, :), allPeaks(ch, :), 'ko'); %ylim([-2 20]);
        title(sprintf('Pearson R %.2f, pVal %.d , Spearman R %.2f, pVal %.d ', ...
            out(ch).allPearsonPeak(1) , out(ch).allPearsonPeak(2),...
            out(ch).allSpearPeak(1) , out(ch).allSpearPeak(2)));
        % save peak stuff
        if allCfg.do_lfpPower2
            figname = sprintf('ch%02d_corrPeak2.png', thisCh);
        else
            figname = sprintf('ch%02d_corrPeak.png', thisCh);
        end
        print(h, fullfile(savefile, figname), '-dpng', '-r300');
        close(h);
    end
elseif strcmp(allCfg.corrLayout, 'overall')
    [rs, ps] = corr([out.stats]', [out.peaks]', 'type', 'Spearman');
    h = figure(); if allCfg.print; set(h, 'visible', 'off'); end;
    plot(allStats, allPeaks, '.');% ylim([-2 20]);
    if strcmp(allCfg.name, 'Hermes')
        title(sprintf('p<0.05 #%d/63, Spearman R %.2f, pVal %.d ',peakCount, rs, ps))
    elseif strcmp(allCfg.name, 'Isis')
        title(sprintf('p<0.05 #%d/22, Spearman R %.2f, pVal %.d ',peakCount, rs, ps))
    end
    figname = sprintf('allChannels_peakCorrelations.png');
    print(h, fullfile(savefile, figname), '-dpng', '-r300');
    close all;
    
    [rs, ps] = corr([out.stats]', [out.rates]', 'type', 'Spearman');
    h = figure(); if allCfg.print; set(h, 'visible', 'off'); end;
    plot(allStats, allRates, '.');% ylim([-2 20]);
    if strcmp(allCfg.name, 'Hermes')
        title(sprintf('p<0.05 #%d/63, Spearman R %.2f, pVal %.d ',rateCount, rs, ps))
    elseif strcmp(allCfg.name, 'Isis')
        title(sprintf('p<0.05 #%d/22, Spearman R %.2f, pVal %.d ',rateCount, rs, ps))
    end
    figname = sprintf('allChannels_rateCorrelations.png');
    print(h, fullfile(savefile, figname), '-dpng', '-r300');
end
