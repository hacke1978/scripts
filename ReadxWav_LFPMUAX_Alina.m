clear all; close all; userpath('clear'); userpath('reset');
% dataDir = '/mnt/hpx/projects/MWNaturalPredict/Hermes/NatImFixMask';
% dataDir = '/mnt/hpx/projects/MWNaturalPredict/isisRaw';
% sesType = 'rfmapping-bar';
% sesType = 'grating-ori';
sesType = 'NatImSEQ';
% sesType = 'NatImFix';
dataDir = '/mnt/archive/MWNaturalPredict/HermesRaw'; %
% dataDir = fullfile('/mnt/hpx/projects/MWNaturalPredict/Hermes', sesType, 'raw');
% dataDir = fullfile('/mnt/v7k/projects/MWNaturalPredict/Hermes', sesType, 'raw');
saveDir = fullfile('/mnt/hpx/projects/MWNaturalPredict/Hermes', sesType);


% saveDir = '/mnt/hpx/projects/MWNaturalPredict/Hermes/grating-ori/';
% saveDir = '/mnt/hpx/projects/MWNaturalPredict/Hermes/NatImFix/';
load('/mnt/hpx/projects/MWNaturalPredict/Cem/matchedChans.mat');
batchTag = '';
% session names
% sessionList = {...
%     %     'hermes_20170531_fixation-naturalim_49',...
%     %     'hermes_20170607_fixation-naturalim_57',...
%     %     'hermes_20170608_fixation-naturalim_58',...
%     %     'hermes_20170608_fixation-naturalim_58b',...
%     %     'hermes_20170612_fixation-naturalim_57',...
%     %     'hermes_20170613_fixation-naturalim_59',...
%     %     'hermes_20170614_fixation-naturalim_61',...
%     %     'hermes_20170616_fixation-naturalim_62',...
%     %     'hermes_20170619_fixation-naturalim_63',...
%     %     'hermes_20170620_fixation-naturalim_64',...
%     %     'hermes_20170620_fixation-naturalim_64b',...
%     %     'hermes_20170623_fixation-naturalim_65',...
%     %     'hermes_20170623_fixation-naturalim_66',...
%     %     'hermes_20170626_fixation-naturalim_66',...
%     %     'hermes_20170626_fixation-naturalim_67',...
%     %     'hermes_20170627_fixation-naturalim_68',...
%     %     'hermes_20170627_fixation-naturalim_68b',...
%     %     'hermes_20170628_fixation-naturalim_69',...
%     %     'hermes_20170628_fixation-naturalim_70',...
%     %     'hermes_20170629_fixation-naturalim_71',...
%     %     'hermes_20170629_fixation-naturalim_72',...
%     %     'hermes_20170630_fixation-naturalim_72',...
%     %     'hermes_20170630_fixation-naturalim_73',...
%     %     'hermes_20170703_fixation-naturalim_73',...
%     %     'hermes_20171120_fixation-naturalim_80',...
%     %     'hermes_20171120_fixation-naturalim_81',...
%     %     'hermes_20171121_fixation-naturalim_82',...
%     %     'hermes_20171121_fixation-naturalim_82b',...
%     %     'hermes_20171122_fixation-naturalim_83',...
%     %     'hermes_20171123_fixation-naturalim_84',...
%     %     'hermes_20171123_fixation-naturalim_84b',...
%     %     'hermes_20171124_fixation-naturalim_85',...
%     %     'hermes_20171127_fixation-naturalim_80',...
%     %     'hermes_20171128_fixation-naturalim_88',...
%     %     'hermes_20171128_fixation-naturalim_89',...
%     %     'hermes_20171129_fixation-naturalim_89',...
%     %     'hermes_20171130_fixation-naturalim_85',...
%     %     'hermes_20171130_fixation-naturalim_89',...
%     %     'hermes_20171201_fixation-naturalim_90',...
%     %     'hermes_20171204_fixation-naturalim_92',...
%     %     'hermes_20171204_fixation-naturalim_93',...
%     %     'hermes_20171205_fixation-naturalim_93',...
%     %     'hermes_20171206_fixation-naturalim_88',...
%     %     'hermes_20171206_fixation-naturalim_89',...
%     %     'hermes_20171207_fixation-naturalim_96',...
%     %     'hermes_20180110_fixation-naturalim_86',...
%     %     'hermes_20180110_fixation-naturalim_96',...
%     %     'hermes_20180118_fixation-naturalim_86B',...
%     %     'hermes_20180118_fixation-naturalim_86G',...
%     %     'hermes_20180118_fixation-naturalim_86R',...
%     %     'hermes_20180122_fixation-naturalim_86B',...
%     %     'hermes_20180122_fixation-naturalim_86G',...
%     %     'hermes_20180122_fixation-naturalim_86k',...
%     %     'hermes_20180122_fixation-naturalim_86R',...
%     %     'hermes_20180124_fixation-naturalim_96b',...
%     %     'hermes_20180123_fixation-naturalim_96',...
%     %     'hermes_20180124_fixation-naturalim_96',...
%     %     'hermes_20180125_fixation-naturalim_97',...
%     %     'hermes_20180125_fixation-naturalim_97b',...
%     %     'hermes_20180125_fixation-naturalim_98',...
%     % 'hermes_20180126_fixation-naturalim_96',...
%     % 'hermes_20180129_fixation-naturalim_101',...
%     % 'hermes_20180129_fixation-naturalim_102',...
%     % 'hermes_20180130_fixation-naturalim_102',...
%     % 'hermes_20180131_fixation-naturalim_102',...
%     'hermes_20180131_fixation-naturalim_103B',...
%     'hermes_20180131_fixation-naturalim_103Y',...
%     };
% sessionList = {
%     'hermes_20170425_fixation-grating-orientation-v2_1',...
%     'hermes_20171208_fixation-grating-orientation-v2_2'
%     };
sessionList = {...
    'hermes_20171113_fixation-color-masksize_78',...
%     'hermes_20171114_fixation-color-masksize_78b',...
%     'hermes_20171115_fixation-color-masksize_79',...
%     'hermes_20171116_fixation-color-masksize_78',...
%     'hermes_20171117_fixation-color-masksize_79b',...
%     'hermes_20170519_fixation-color-masksize_42'
    };
% sessionList = {
%     'hermes_20180118_fixation-naturalim_86B',...
%     'hermes_20180118_fixation-naturalim_86G',...
%     'hermes_20180118_fixation-naturalim_86R',...
%     'hermes_20180122_fixation-naturalim_86B',...
%     'hermes_20180122_fixation-naturalim_86G',...
%     'hermes_20180122_fixation-naturalim_86k',...
%     'hermes_20180122_fixation-naturalim_86R',...
% };
% sessionList = {
% 'hermes_20180201_fixation-naturalim_103E',...
% 'hermes_20180201_fixation-naturalim_103K',...
% 'hermes_20180202_fixation-naturalim_103E',...
% 'hermes_20180202_fixation-naturalim_103G',...
% 'hermes_20180202_fixation-naturalim_103R',...
% };
% sessionList = {
%     'hermes_20170703_fixation-naturalim_74',...
%     'hermes_20170704_fixation-naturalim_75',...
%     'hermes_20170706_fixation-naturalim_76',...
%     'hermes_20170706_fixation-naturalim_77',...
%     'hermes_20170712_fixation-naturalim_60',...
%     'hermes_20170724_fixation-naturalim_70',...
%     'hermes_20170725_fixation-naturalim_71',...
%     };
% sessionList = {
% 'hermes_20180205_fixation-naturalim_102',...
% 'hermes_20180207_fixation-naturalim_103E',...
% 'hermes_20180207_fixation-naturalim_103G3',...
% 'hermes_20180207_fixation-naturalim_103K',...
% 'hermes_20180207_fixation-naturalim_103W'
% };
% sessionList = {
%     'hermes_20180208_fixation-naturalim_103B2',...
%     'hermes_20180208_fixation-naturalim_103R2',...
%     'hermes_20180208_fixation-naturalim_103Y2',...
%     'hermes_20180209_fixation-naturalim_103G4',...
%     'hermes_20180209_fixation-naturalim_103W2',...
%     'hermes_20180213_fixation-naturalim_104',...
%     'hermes_20180213_fixation-naturalim_105',...
%     'hermes_20180213_fixation-naturalim_106',...
%     };
% sessionList = {
% 'hermes_20170522_fixation-color-masksize_43b',...
% };
% sessionList = {
%     %     'hermes_20180129_fixation-naturalim_101',...
% %     'hermes_20180129_fixation-naturalim_102',...
% %     'hermes_201800129_fixation-naturalim_101',...
% %     'hermes_20180126_fixation-naturalim_96'
%     };
%% Functions
failedSessions = [];
for ses = 1:length(sessionList)
    %     try
    global outputDirSes
    outputDirSes = fullfile(saveDir, strcat(sessionList{ses}));
    mkdir(outputDirSes); %mkdir(fullfile(outputDirSes, 'data'));
    
    if strcmp(sessionList{ses}, 'hermes_20171113_fixation-color-masksize_78')
        dataDirSes = fullfile(dataDir, 'hermes_20171113_fixation-naturalim-masked_78');
    else
        dataDirSes = fullfile(dataDir, sessionList{ses});
    end
    
    
    % Load neural data - LFP, MUA, Waveforms
    if strcmp(sessionList{ses}, 'hermes_20171113_fixation-color-masksize_78')
        pathXWav = dir(fullfile(dataDirSes, sprintf('%s_xWav*', 'hermes_20171113_fixation-naturalim-masked_78')));
    else
        pathXWav = dir(fullfile(dataDirSes, sprintf('%s_xWav*', sessionList{ses})));
    end
    sesName = sessionList{ses};
    %% RAW to MUAX / LFP
    %     pathLfp = dir(fullfile(outputDirSes, '*_xWav.lfp'));
    %     pathMuax = dir(fullfile(outputDirSes, '*_xWav.muax'));
    %     pathMuat = dir(fullfile(outputDirSes, '*_xWav.muat'));
    %         if isempty(pathLfp) || isempty(pathMuax)
%     cfg = [];
%     cfg.filename = strcat(dataDirSes, '/', {pathXWav.name});
%     cfg.targetFolder = outputDirSes;
%     cfg.calcLocation = 'slurm';
%     tdt_preprocessing_AP(cfg);
%     %     % %         %     end
%     %     % % %         if isempty(pathMuat)
%     cfg = [];
%     cfg.filename = strcat(dataDirSes, '/', {pathXWav.name});
%     cfg.targetFolder = outputDirSes;
%     cfg.calcLocation = 'slurm';
%     cfg.pSigma = 3;
%     cfg.pISI = 1.5;
%     %     cfg.pTime = [0 60]; % in sec
%     tdt_extractspikes_AP(cfg);
    %     %             end
    %% Chop into Trials
    pathLfp = dir(fullfile(outputDirSes, '*_xWav.lfp'));
    pathMua = dir(fullfile(outputDirSes, '*_xWav.mua'));
    pathMuat = dir(fullfile(outputDirSes, '*_xWav.muat'));
    pathMuax = dir(fullfile(outputDirSes, '*On.muax'));
    %     pathMuax = dir(fullfile(outputDirSes, '*_xWav.muax'));
    load(fullfile(outputDirSes, pathLfp.name), '-mat');
    lfp.data = data; clear data;
    load(fullfile(outputDirSes, pathMuax.name), '-mat');
    muax.data = data; clear data;
    load(fullfile(outputDirSes, pathMuat.name), '-mat');
    muat.data = data; clear data;
    load(fullfile(outputDirSes, pathMua.name), '-mat');
    mua.data = data; clear data;
    %% load TDT strobes
    if strcmp(sesType, 'grating-ori')
        preStim = 1.0;
    elseif strcmp(sesType, 'rfmapping-bar')
        preStim = 0;
    else
        preStim = 0.5;
    end
    cfg = [];
    cfg.trl = muax.data.cfg.trl;
    trialinfo = muax.data.trialinfo;
    sampleinfo = muax.data.sampleinfo;
    label_tdt = muax.data.label_tdt(matchedChans);
    label = muax.data.label; % get labels
    [~, channelSortInd] = sort(matchedChans);
    
    % get stuff
    Cond = muax.data.trialinfo(:, 2);
    trialStart = muax.data.sampleinfo(:, end-1)/muax.data.fsample;
    trialEnd = muax.data.sampleinfo(:, end)/muax.data.fsample;
    sesName = sessionList{ses};
    saveName = sprintf('/mnt/hpx/projects/MWNaturalPredict/Hermes/%s/%s/%s', sesType, sesName, fullfile(sesName, batchTag));
    
    if strcmp(sesType, 'NatImSEQ')
        maskFirst = trialinfo(:, end-5);
        maskSize = trialinfo(:, end-4);
        maskPos = trialinfo(:, end-2);
        maskSize(maskSize==6 & maskPos==2) = 6.01;
        maskSize(maskSize==6.01) = 4; maskSize(maskSize==2.0) = 3;
        maskSize(maskSize==1.0) = 2; maskSize(maskSize==0.5) = 1;
        maskLib = unique(maskSize);
        seqLib = unique(maskFirst);
        condLib = unique(Cond);
        Cond = maskFirst*(length(condLib)*length(maskLib))+length(maskLib)*(Cond-1)+maskSize;
    end
    
    if strcmp(sesType, 'rfmapping-bar')
        taccept = true(size(muax.data.trialinfo, 1), 1);
        trialinfo = [[1:length(Cond)]' taccept Cond trialStart trialEnd];
    elseif strcmp(sesType, 'grating-ori')
        taccept = ((muax.data.trialinfo(:, 4)));
        CondSF = muax.data.trialinfo(:, 3);
        trialinfo = [[1:length(Cond)]' taccept Cond CondSF trialStart trialEnd];
    else
        taccept = ((muax.data.trialinfo(:, 8)));
        trialinfo = [[1:length(Cond)]' taccept Cond trialStart trialEnd];
    end
    
    
    pathMuax = dir(fullfile(outputDirSes, '*_xWav.muax'));
    load(fullfile(outputDirSes, pathMuax.name), '-mat');
    muax.data = data; clear data;
    
    %% Chop it - LFP
    data.fsample = lfp.data.fsample;
    data.cfg = cfg;
    data.sampleinfo = sampleinfo;
    %     data.cfg.trl = muax.data.cfg.trl;
    data.label_tdt = lfp.data.label(matchedChans);
    data.label = label;
    data.trialinfo = trialinfo;% add in trialno
    %         data.trialinfo = trialinfo;
    time = lfp.data.time{1};
    for ii = 1:length(Cond)
        %         trialChosen = trialStart(ii) <= time & time <= trialEnd(ii);
        %         data.trial{ii} = lfp.data.trial{1}(matchedChans, trialChosen);
        data.trial{ii} = lfp.data.trial{1}(matchedChans, sampleinfo(ii, 1):sampleinfo(ii, 2));
        %         data.time{ii} = time(trialChosen)-trialStart(ii)-preStim;
        data.time{ii} = time(sampleinfo(ii, 1):sampleinfo(ii, 2))-time(sampleinfo(ii, 1))-preStim;
    end
    trialtime = [cellfun(@min, data.time)' cellfun(@max, data.time)'];
    save([saveName '_chopped.lfp'],'data', '-v7.3');
    clear data
    
    %% Chop it - MUAX
    %     preStim = 0.5;
    data.fsample = muax.data.fsample;
    data.cfg = cfg;
    data.sampleinfo = sampleinfo;
    %     data.cfg.trl = muax.data.cfg.trl;
    %     data.label_tdt = muax.data.label;
    data.label_tdt = muax.data.label(matchedChans);
    data.label = label;
    data.trialinfo = trialinfo;
    %         data.trialinfo = trialinfo;
    time = muax.data.time{1};
    for ii = 1:length(Cond)
        trialChosen = trialStart(ii) <= time & time <= trialEnd(ii);
        %         data.trial{ii} = muax.data.trial{1}(matchedChans, trialChosen);
        data.trial{ii} = muax.data.trial{1}(matchedChans, sampleinfo(ii, 1):sampleinfo(ii, 2));
        %         data.time{ii} = time(trialChosen)-trialStart(ii)-preStim;
        data.time{ii} = time(sampleinfo(ii, 1):sampleinfo(ii, 2))-time(sampleinfo(ii, 1))-preStim;
    end
    save([saveName '_chopped.muax'],'data', '-v7.3');
    clear data
    
    %% Chop it - MUA
    data.fsample = mua.data.fsample;
    data.cfg = cfg;
    data.sampleinfo = sampleinfo;
    data.label_tdt = mua.data.label(matchedChans);
    data.label = label;
    data.trialinfo = trialinfo;
    %         data.trialinfo = trialinfo;
    time = mua.data.time{1};
    for ii = 1:length(Cond)
        trialChosen = trialStart(ii) <= time & time <= trialEnd(ii);
        %         data.trial{ii} = muax.data.trial{1}(matchedChans, trialChosen);
        data.trial{ii} = mua.data.trial{1}(matchedChans, sampleinfo(ii, 1):sampleinfo(ii, 2));
        %         data.time{ii} = time(trialChosen)-trialStart(ii)-preStim;
        data.time{ii} = time(sampleinfo(ii, 1):sampleinfo(ii, 2))-time(sampleinfo(ii, 1))-preStim;
    end
    save([saveName '_chopped.mua'],'data', '-v7.3');
    clear data
    
    %% Chop it - MUAThreshold
    %     preStim = 0.5;
    data.fsample = muat.data.fsample;
    data.cfg = muat.data.cfg;
    data.sampleinfo = muat.data.sampleinfo;
    data.label = muat.data.label;
    %     data.trialinfo = muax.data.trialinfo;
    data.trialinfo = trialinfo;
    %         data.trialinfo = trialinfo;
    time = [data.sampleinfo(1):data.sampleinfo(2)]/data.fsample;
    
    % Make Spikes
    spike.hdr = muat.data.hdr;
    %         spike.trialinfo = data.trialinfo;
    %     spike.sampleinfo = muax.data.trialinfo;
    spike.trialinfo = trialinfo;
    nrTrials = length(Cond);
    for ch = 1:length(muat.data.label)
        spikeTimes = muat.data.trial{ch};
        spikeStamps = round(data.fsample*muat.data.trial{ch})-1;
        thisTDTCh = find(matchedChans == (str2num(muat.data.label{ch}(6:end))));
        thisIdx = find(strncmp(label_tdt(channelSortInd), sprintf('aMUA%d', thisTDTCh), 10));
        data.sampleinfo = [];
        for ii = 1:nrTrials
            data.sampleinfo = [data.sampleinfo; round(data.fsample*[trialStart(ii) trialEnd(ii)])];
            data.trial{thisIdx, ii} = spikeTimes(and(spikeTimes>trialStart(ii), spikeTimes<=trialEnd(ii)))- trialStart(ii) - preStim;
            data.timestamp{thisIdx, ii} = spikeStamps(and(spikeTimes>trialStart(ii), spikeTimes<=trialEnd(ii)));
        end
        spike.timestamp{thisIdx} = ([data.timestamp{thisIdx, :}]');
        spike.time{thisIdx} = ([data.trial{thisIdx, :}]');
        % fix the channel order + labels to match muax & lfp
        spike.label{thisIdx} = label{thisIdx};
        spike.label_tdt{thisIdx} = muat.data.label{ch};
        % fix the channel order + labels to match muax & lfp
        data.label{thisIdx} = label{thisIdx};
        data.label_tdt{thisIdx} = muat.data.label{ch};
        % done fixing
        spike.sampleinfo = data.sampleinfo;
        spike.trialtime = trialtime;
        %                 spike.trialtime = repmat([-0.5 1.3], nrTrials, 1);
        %         spike.trialtime = repmat([-0.5 3.8], nrTrials, 1);
        %         spike.trialtime = repmat([0 3.8], nrTrials, 1);
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
    
    save([saveName '_chopped.muat'],'data', '-v7.3');
    clear data
    
    save([saveName '_chopped.spike'],'spike', '-v7.3');
    clear spike
    % % %
    % %     %     %% Get the eye data
    % %     %     if str2num(sesName(6:8)) > 458
    % %     %         % Convert the EDF file to ASC
    % %     %         pathEDF = dir(fullfile(outputDirSes, '*.edf'));
    % %     %         system(sprintf('linuxinstall/bin/edf2asc -t -c -sp %s', fullfile(outputDirSes, pathEDF.name)));
    % %     %
    % %     %         % Get the ASC
    % %     %         eyeList = dir(fullfile(outputDirSes, '*Eye.asc'));
    % %     %         cfg = [];
    % %     %         cfg.dataset = fullfile(outputDirSes, eyeList.name);
    % %     %         cfg.trialOutcome = taccept;
    % %     %         cfg.savepath = saveName;
    % %     %         dataEye = convert_edf_fieldtrip(cfg, cfg.dataset);
    % %     %
    % %     %         if sum(dataEye.taccept ~= cfg.trialOutcome)
    % %     %             error('Mismatch conditions eye data')
    % %     %         end
    % %     %         %% Make trial selection
    % %     %         preStim = 0.5;
    % %     %         data.hdr = dataEye.hdr;
    % %     %         data.label = dataEye.label;
    % %     %         data.fsample = dataEye.fsample;
    % %     %         data.sampleinfo= round(dataEye.fsample*[trialStart trialEnd]); % get sample number
    % %     %         data.trialinfo = [[1:length(Cond)]' taccept' Cond trialStart trialEnd]; % add in trialno
    % %     %         %     time  = 0:1/data.fsample:(1/data.fsample)*(size(Eyes, 2)-1);
    % %     %         time = dataEye.time{1};
    % %     %         for ii = 1:length(Cond)
    % %     %             trialChosen = trialStart(ii) <= time & time <= trialEnd(ii);
    % %     %             data.trial{ii} = dataEye.trial{1}(:, trialChosen);
    % %     %             data.time{ii} = time(trialChosen)-trialStart(ii)-preStim;
    % %     %         end4
    % %     %         save([MATFilename '.eyes'],'data');
    % %     %         clear data
    % %     %     end
    clearvars -except dataDir saveDir batchTag sessionList sesType matchedChans failedSessions sesName
    %     catch
    %         failedSessions = [failedSessions; {sesName}];
    %     end
end
