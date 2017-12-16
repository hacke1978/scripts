clear all; close all; userpath('clear'); userpath('reset');
dataDir = '/mnt/archive/MWNaturalPredict';
saveDir = '/mnt/hpx/projects/MWNaturalPredict/Hermes/NatImFix/';
load('/mnt/hpx/projects/MWNaturalPredict/Cem/matchedChans.mat');
batchTag = '';

% session names
sessionDir = {
'hermes_20170808_rfmapping-bar_1'
    };
    %% Functions
    for ses = 1:length(sessionDir)
    global outputDirSes
    outputDirSes = fullfile(saveDir, strcat(batchTag, sessionDir{ses}));
    mkdir(outputDirSes); mkdir(fullfile(outputDirSes, 'data'));
    dataDirSes = fullfile(dataDir, sessionDir{ses});
    % Load neural data - LFP, MUA, Waveforms
    pathXWav = dir(fullfile(dataDirSes, '*xWav*'));
    
    %% RAW to MUAX / LFP
    cfg = [];
    cfg.filename = strcat(dataDirSes, '/', {pathXWav.name});
    cfg.targetFolder = outputDirSes;
    cfg.calcLocation = 'slurm';
    tdt_extractspikes_AP(cfg);
    
    %% Chop into Trials
    pathLfp = dir(fullfile(outputDirSes, '*.lfp'));
    pathMuax = dir(fullfile(outputDirSes, '*.muax'));
    pathMuat = dir(fullfile(outputDirSes, '*_xWav.muat'));
    load(fullfile(outputDirSes, pathLfp.name), '-mat');
    lfp.data = data; clear data;
    load(fullfile(outputDirSes, pathMuax.name), '-mat');
    muax.data = data; clear data;
    load(fullfile(outputDirSes, pathMuat.name), '-mat');
    muat.data = data; clear data;
    
    %% load TDT strobes
    Cond = muax.data.trialinfo(:, 2);
    taccept = ones(size(Cond)); % all correct trials
    trialStart = muax.data.sampleinfo(:, end-1)/muax.data.fsample;
    trialEnd = muax.data.sampleinfo(:, end)/muax.data.fsample;
    sesName = sessionDir{ses};
    saveName = fullfile(saveDir, sprintf('%s/%s', sesName, sesName));
  
    % Chop it - MUAThreshold
    preStim = 0.0;
    data.fsample = muat.data.fsample;
    data.cfg = muat.data.cfg;
    data.sampleinfo = muat.data.sampleinfo;
    data.label = muat.data.label;
    %     data.trialinfo = muax.data.trialinfo;
    data.trialinfo = [[1:length(Cond)]' taccept Cond trialStart trialEnd]; % add in trialno
    time = [data.sampleinfo(1):data.sampleinfo(2)]/data.fsample;
    
    % Make Spikes
    spike.hdr = muat.data.hdr;
    spike.trialinfo = muax.data.trialinfo;
    %     spike.sampleinfo = muax.data.trialinfo;
    %     spike.trialinfo = [[1:length(Cond)]' Cond trialStart trialEnd]; % add in trialno
    nrTrials = length(Cond);
    for ch = 1:length(muat.data.label)
        spikeTimes = muat.data.trial{ch};
        thisTDTCh = find(matchedChans == (str2num(muat.data.label{ch}(6:end))));
        thisIdx = find(strncmp(muax.data.label_tdt, sprintf('aMUA%d', thisTDTCh), 10));
        data.sampleinfo = [];
        for ii = 1:nrTrials
            data.sampleinfo = [data.sampleinfo; round(data.fsample*[trialStart(ii) trialEnd(ii)])];
            data.trial{thisIdx, ii} = spikeTimes(and(spikeTimes>trialStart(ii), spikeTimes<=trialEnd(ii)))- trialStart(ii) - preStim;
        end
        spike.time{thisIdx} = ([data.trial{thisIdx, :}]');
        % fix the channel order + labels to match muax & lfp
        spike.label{thisIdx} = muax.data.label{thisIdx};
        spike.label_tdt{thisIdx} = muat.data.label{ch};
        % fix the channel order + labels to match muax & lfp
        data.label{thisIdx} = muax.data.label{thisIdx};
        data.label_tdt{thisIdx} = muat.data.label{ch};
        % done fixing
        spike.sampleinfo = data.sampleinfo;
        spike.trialtime = repmat([-0.5 1.2], nrTrials, 1);
        spike.cond = Cond;
        spike.taccept = taccept;
        onesVec = cellfun(@(x) (ones(length(x), 1)'), data.trial(thisIdx, :), 'UniformOutput', false);
        trialVec = cellfun(@(x, y) (x*y), onesVec, num2cell([1:nrTrials]), 'UniformOutput', false);
        spike.trial{thisIdx} = [trialVec{:}]';
    end
    spike.sampleinfo = data.sampleinfo;
    spike.label = spike.label';
    spike.label_tdt = spike.label_tdt';
    data.label_tdt = data.label_tdt';
    
    save([saveName '_chopped.muat'],'data')
    clear data
    
    save([saveName '_chopped.spike'],'spike')
    clear spike
    

    clearvars -except dataDir saveDir batchTag sessionDir
    
    end
