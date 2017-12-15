function [out] =  getRatesAndStats(allCfg, allFiles)
% this script correlates the rates and the im stats
out = [];
statDir = '/mnt/hpx/projects/MWNaturalPredict/Cem/Analysis/Inpainting';
stimSetAll = '/mnt/v7k/home/uranc/workspace/VinckLab/Dataset/Processed/_6flickr/stimsetChosen';
load(fullfile(stimSetAll, 'imSetHighSinger.mat'));
load(fullfile(stimSetAll, 'imSetLowSinger.mat'));
allStats = horzcat(highStats, lowStats);
if strcmp(allCfg.name, 'hermes')
    load(fullfile(statDir, 'allChSSIM_hermes.mat'));
%     load(fullfile('/mnt/hpx/projects/MWNaturalPredict/Cem/Analysis/RFmaps/', sprintf('GRFs_Fit_%02d.mat', 64)));
    RFs = allFiles(1).RFs;
    label = [RFs.label];
    label = cellfun(@(x) (x(4:end)), label, 'UniformOutput', false);
    caccept = ones(1, 64);
    if allCfg.isOri
        % get channels ori
        oriRFs = [RFs.ori];
        
        % get im ori
        load(fullfile(statDir, sprintf('allOri_%s.mat', allCfg.name)))
    end
elseif strcmp(allCfg.name, 'isis')
    load(fullfile('/mnt/hpx/projects/MWNaturalPredict/Cem/Analysis/RFmaps/', sprintf('GRFs_Fit_v1.mat')))
    load(fullfile(statDir, 'allChSSIM_isis.mat'));
    caccept = ~ismember(1:32, [9 13 14 16 22 23 24 26 31 32]);
    label = strsplit(num2str(find(caccept)));
end

if isfield(allFiles, 'trialPSTH') && isfield(allFiles, 'lfpPower')
    
    % first get natural images
    allNames = vertcat({allChSSIM(1, :).stimName});
    stimList = vertcat({allFiles.sname});
    imInd = cellfun(@(x) find(strncmp(allNames, x, 15)), stimList, 'UniformOutput', false);
    grayInd = cellfun(@(x) isempty(x), imInd, 'UniformOutput', false);
    condSelect = ~[grayInd{:}]; % get natural images only
    allNames = vertcat({allChSSIM(1, :).stimName});
    stimList = vertcat({allFiles.sname});
    imInd = cellfun(@(x) find(strncmp(allNames, x, 15)), stimList, 'UniformOutput', false);
    imInd = [imInd{:}];
    
    % get all psth
    allPSTH = cat(3, allFiles(condSelect).trialPSTH);
    ttime = allPSTH(1).time;
    tTransient = allCfg.window_rate(1, 1) < ttime & ttime < allCfg.window_rate(1, 2);
    tSustained = allCfg.window_rate(2, 1) < ttime & ttime < allCfg.window_rate(2, 2);
    allThis = cat(3, allPSTH.avg);
    
    % get all power
    allPowerLFP = cat(3, allFiles(condSelect).lfpPower);
    tfreq = allPowerLFP(1).freq;
    tFreqLow = allCfg.window_freq(1, 1) < tfreq & tfreq < allCfg.window_freq(1, 2);
    tFreqHigh = allCfg.window_freq(2, 1) < tfreq & tfreq < allCfg.window_freq(2, 2);
    data = cat(3, allPowerLFP.powspctrm);
    if isfield(allFiles, 'lfpPower2')
        powerLFP2 = cat(2, allFiles(condSelect).lfpPower2);
        data2 = cat(3, powerLFP2.powspctrm);
        data = (data + data2)/2;
        disp('averaging two windows!!!');
    end
    
    % get baseline
    powerBaseline = cat(2, allFiles(1).powerBaseline);
    baseline = cat(3, powerBaseline.powspctrm);
    
    % ready for the correlations
    allSpearRate = []; allPearsonRate = [];
    allSpearPeak = []; allPearsonPeak = [];
    for ch=1:size(allChSSIM, 1)
        % fsim values
        predValues = [allChSSIM(ch, imInd).fsim];
        
        % now get the rates
        allTransient = mean(squeeze(allThis(ch, tTransient, :)), 1);
        allSustained = mean(squeeze(allThis(ch, tSustained, :)), 1);
        allRates = abs(allTransient./(1e-20 + abs(allSustained)));
        [sCorr, sVal] = corr(predValues', allRates', 'type', 'Spearman');
        %         allSpearRate = [allSpearRate; [sCorr, sVal]];
        allSpearRate = [sCorr, sVal];
        [pCorr, pVal] = corr(predValues', allRates', 'type', 'Pearson');
        %         allPearsonRate = [allPearsonRate; [pCorr, pVal]];
        allPearsonRate = [pCorr, pVal];
        % now get the gamma
        % normalized to baseline
        powerToBase = (log10(squeeze(data(ch, :, :))./squeeze(repmat(1e-20+baseline(ch, :), 1, 1, sum(condSelect)))));
        
        allLow = max(powerToBase(tFreqLow, :), [], 1);
        allHigh = mean(powerToBase(tFreqHigh, :), 1);
        allPeak = allLow./(1e-20 + abs(allHigh));
        [sCorr, sVal] = corr(predValues', allPeak', 'type', 'Spearman');
        %         allSpearPeak = [allSpearPeak; [sCorr, sVal]];
        allSpearPeak = [sCorr, sVal];
        [pCorr, pVal] = corr(predValues', allPeak', 'type', 'Pearson');
        %         allPearsonPeak = [allPearsonPeak; [pCorr, pVal]];
        allPearsonPeak = [pCorr, pVal];
        
        
%           if allCfg.isOri
%             oriCond = [allOri(ch, imInd).peakOri];
%             deltaOri = abs(oriCond-repmat(oriRFs(ch), 1, length(oriCond)));
%             deltaOri(deltaOri>90) = deltaOri(deltaOri>90)-90;
%             deltaOri = deltaOri/90;
            
% % %             [a, b] = corr(deltaOri', allPeak')
% %             X = [predValues; deltaOri]';
% % %             X = [predValues]';
% %             y = allPeak;
% % %             y = allRates;
% %             
% %             [bPred, devPred, statsPred] = glmfit(predValues, y);
% %             [bOri, devOri, statsOri]  = glmfit(deltaOri, y);
% %             [bFull, devFull, statsFull]  = glmfit(X, y);
% %             yhatPred = glmval(bPred, X, 'identity');
% %             yhatOri = glmval(bOri, X, 'identity');
% %             yhatFull = glmval(bFull, X, 'identity');
% % 
% % 
% %             fitlm(X, y)
% %             
% %             figure
% %             plot(X, y, '.'); hold on
% %             plot(X, yhat, 'r.')
%         end
        
        
        % outs
        out(ch).rates = allRates;
        out(ch).peaks = allPeak;
        out(ch).stats = predValues;
        out(ch).allSpearRate = allSpearRate;
        out(ch).allPearsonRate = allPearsonRate;
        out(ch).allSpearPeak = allSpearPeak;
        out(ch).allPearsonPeak = allPearsonPeak;
        out(ch).label = label(ch);
        if allCfg.isOri
        out(ch).ori = deltaOri;
        end
    end
    
elseif isfield(allFiles, 'timelockMUAX')
else
    warning('no rate files!')
end
end