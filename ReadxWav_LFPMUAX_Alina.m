clear all; close all; userpath('clear'); userpath('reset');
% dataDir = '/mnt/hpx/projects/MWNaturalPredict/Hermes/NatImFixMask';
% dataDir = '/mnt/hpx/projects/MWNaturalPredict/isisRaw';
dataDir = '/mnt/archive/MWNaturalPredict';
saveDir = '/mnt/hpx/projects/MWNaturalPredict/Hermes/NatImFix/';
load('/mnt/hpx/projects/MWNaturalPredict/Cem/matchedChans.mat');
batchTag = '5sigma';
% session names
sessionDir = {
    ...hermes_20171120_fixation-naturalim_80/%'hermes_20170524_fixation-naturalim_13', ...
    ...%'hermes_20170529_fixation-naturalim_16', ...
    ...%'hermes_20170529_fixation-naturalim_17',...
    ...% 'hermes_20171121_fixation-naturalim_82'
    ...%     'hermes_20171121_fixation-naturalim_82',...
%     'hermes_20171123_fixation-naturalim_84',...
%     'hermes_20171120_fixation-naturalim_80',...
%     'hermes_20171124_fixation-naturalim_85',...
%     'hermes_20171127_fixation-naturalim_86',...
%     'hermes_20171127_fixation-naturalim_87',...
%     'hermes_20171128_fixation-naturalim_88',...
%     'hermes_20171201_fixation-naturalim_90',...
%     'hermes_20171204_fixation-naturalim_92'
% 'hermes_20171121_fixation-naturalim_82',...
% 'hermes_20171122_fixation-naturalim_83'
% 'hermes_20171130_fixation-naturalim_89',...
% 'hermes_20171207_fixation-naturalim_96',...
% 'hermes_20171205_fixation-naturalim_93',...
    'hermes_20170808_rfmapping-bar_1'
    };
% sessionDir = dir(fullfile(dataDir, '*fixation-naturalim*'));
% sessionDir = dir(fullfile(dataDir, '*fixation-grating-ori*'));
%% Functions
for ses = 1:length(sessionDir)
    global outputDirSes
    outputDirSes = fullfile(saveDir, strcat(sessionDir{ses}));
    mkdir(outputDirSes); mkdir(fullfile(outputDirSes, 'data'));
    dataDirSes = fullfile(dataDir, sessionDir{ses});
    % Load neural data - LFP, MUA, Waveforms
    pathXWav = dir(fullfile(dataDirSes, '*xWav*'));
    
    %% RAW to MUAX / LFP
%     cfg = [];
%     cfg.filename = strcat(dataDirSes, '/', {pathXWav.name});
%     cfg.targetFolder = outputDirSes;
%     cfg.calcLocation = 'slurm';
    %     tdt_preprocessing_AP(cfg);
    
    cfg = [];
    cfg.filename = strcat(dataDirSes, '/', {pathXWav.name});
    cfg.targetFolder = outputDirSes;
    cfg.calcLocation = 'slurm';
    cfg.pSigma = 5;
    cfg.pISI = 1.5;
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
%         taccept = muax.data.trialinfo(:, 8);
    Cond = muax.data.trialinfo(:, 2);
    taccept = ones(size(Cond));
    trialStart = muax.data.sampleinfo(:, end-1)/muax.data.fsample;
    trialEnd = muax.data.sampleinfo(:, end)/muax.data.fsample;
    sesName = sessionDir{ses};
    saveName = sprintf('/mnt/hpx/projects/MWNaturalPredict/Hermes/NatImFix/%s/%s', sesName, fullfile(sesName, batchTag));
%         [trialEnd, trialStart, Cond, taccept] = getSnipsForSession(sesName, 'NatIm');
    
    % %     % Chop it - LFP
    % %     preStim = 0.5;
    % %     data.fsample = lfp.data.fsample;
    % %     data.cfg = lfp.data.cfg;
    % %     data.sampleinfo = lfp.data.sampleinfo;
    % %     data.label = lfp.data.label;
    % %     data.trialinfo = [[1:length(Cond)]' taccept Cond trialStart trialEnd]; % add in trialno
    % %     time = lfp.data.time{1};
    % %     for ii = 1:length(Cond)
    % %         trialChosen = trialStart(ii) <= time & time <= trialEnd(ii);
    % %         data.trial{ii} = lfp.data.trial{1}(:, trialChosen);
    % %         data.time{ii} = time(trialChosen)-trialStart(ii)-preStim;
    % %     end
    % %     save([saveName '_chopped.lfp'],'data');
    % %     clear data
    % %
    % %     % Chop it - MUAX
    % %     preStim = 0.5;
    % %     data.fsample = muax.data.fsample;
    % %     data.cfg = muax.data.cfg;
    % %     data.sampleinfo = muax.data.sampleinfo;
    % %     data.label = muax.data.label;
    % %     data.trialinfo = [[1:length(Cond)]' taccept Cond trialStart trialEnd]; % add in trialno
    % %     time = muax.data.time{1};
    % %     for ii = 1:length(Cond)
    % %         trialChosen = trialStart(ii) <= time & time <= trialEnd(ii);
    % %         data.trial{ii} = muax.data.trial{1}(:, trialChosen);
    % %         data.time{ii} = time(trialChosen)-trialStart(ii)-preStim;
    % %     end
    % %     save([saveName '_chopped.muax'],'data');
    % %     clear data
    
    % Chop it - MUAThreshold
    preStim = 0.5;
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
    
    %     %% Get the eye data
    %     if str2num(sesName(6:8)) > 458
    %         % Convert the EDF file to ASC
    %         pathEDF = dir(fullfile(outputDirSes, '*.edf'));
    %         system(sprintf('linuxinstall/bin/edf2asc -t -c -sp %s', fullfile(outputDirSes, pathEDF.name)));
    %
    %         % Get the ASC
    %         eyeList = dir(fullfile(outputDirSes, '*Eye.asc'));
    %         cfg = [];
    %         cfg.dataset = fullfile(outputDirSes, eyeList.name);
    %         cfg.trialOutcome = taccept;
    %         cfg.savepath = saveName;
    %         dataEye = convert_edf_fieldtrip(cfg, cfg.dataset);
    %
    %         if sum(dataEye.taccept ~= cfg.trialOutcome)
    %             error('Mismatch conditions eye data')
    %         end
    %         %% Make trial selection
    %         preStim = 0.5;
    %         data.hdr = dataEye.hdr;
    %         data.label = dataEye.label;
    %         data.fsample = dataEye.fsample;
    %         data.sampleinfo= round(dataEye.fsample*[trialStart trialEnd]); % get sample number
    %         data.trialinfo = [[1:length(Cond)]' taccept' Cond trialStart trialEnd]; % add in trialno
    %         %     time  = 0:1/data.fsample:(1/data.fsample)*(size(Eyes, 2)-1);
    %         time = dataEye.time{1};
    %         for ii = 1:length(Cond)
    %             trialChosen = trialStart(ii) <= time & time <= trialEnd(ii);
    %             data.trial{ii} = dataEye.trial{1}(:, trialChosen);
    %             data.time{ii} = time(trialChosen)-trialStart(ii)-preStim;
    %         end
    %         save([MATFilename '.eyes'],'data');
    %         clear data
    %     end
    clearvars -except dataDir saveDir batchTag sessionDir matchedChans
    
end
