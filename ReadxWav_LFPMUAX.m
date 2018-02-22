clear all; close all; userpath('clear');
addpath('/opt/ESIsoftware/slurmfun/')
addpath('/mnt/hpx/opt/ESIsoftware/matlab/')
addpath('/mnt/hpx/slurm/uranc/fieldtrip/')
ft_defaults
%% Dirs
sesType = 'NatImFix'; %'NatImSEQ' 'NatImFix' 'rfmapping-bar'
batchTag = '';

% session names
sessionList = {
%     'ares025a03',...
%     'ares026a01',...
%     'ares026a02',...
%     'ares027a01',...
%     'ares027a02',...
%     'ares028a01',...
%     'ares029a01',...
%     'ares030a01',...
%     'ares030a02',...
%     'ares033a01',...
%     'ares033a02',...
%     'ares033a03',...
%     'ares034a03',...
%     'ares034a04',...
%     'ares037a01',...
%     'ares038a01',...
%     'ares038a02',...
%     'ares039a01',...
%     'ares039a02',...
% 'ares040a01',...
'ares040a02',...
% 'ares040a03',...
% 'ares042a01',...
% 'ares042a02',...
    };
% sessionList = {
%     'ares023a02',...
%     'ares024a01',...
%     'ares025a01',...
%     'ares025a02',...
%     };
% sessionList = {'ares027a02'};
%% Functions

failedSessions = [];
for ses = 1:length(sessionList)
    global outputDirSes
    %     try
    % get monkey dirs
    monkeyName = sessionList{ses}(1:4);
    sesName = sessionList{ses};
    if strcmp(monkeyName, 'herm'); monkeyName = 'Hermes'; end
    if strcmp(monkeyName, 'ares'); monkeyName = 'Ares'; end
    if strcmp(monkeyName, 'isis'); monkeyName = 'Isis'; end
    
    % dir
    dataDir = fullfile('/mnt/archive/MWNaturalPredict', strcat(monkeyName, 'Raw'));
    ESIsaveDir = fullfile('/mnt/hpx/projects/MWNaturalPredict', monkeyName, sesType);
    outputDirSes = fullfile(ESIsaveDir, strcat(batchTag, sessionList{ses}));
    mkdir(outputDirSes); mkdir(fullfile(outputDirSes, 'data'));
    dataDirSes = fullfile(dataDir, sessionList{ses});
    pathXWav = dir(fullfile(dataDirSes, '*.sev'));  % get raw data
    %% RAW to MUAX / LFP
%     cfg = [];
%     cfg.filename = strcat(dataDirSes, '/', {pathXWav.name});
%     cfg.targetFolder = outputDirSes;
%     cfg.calcLocation = 'slurm';
%     tdt_preprocessing_AP(cfg);
%     %
%     cfg = [];
%     cfg.filename = strcat(dataDirSes, '/', {pathXWav.name});
%     cfg.targetFolder = outputDirSes;
%     cfg.calcLocation = 'slurm';
%     cfg.pSigma = 3; % amp. threshold pSigma*median(abs(x)/0.6745)) % Quian Quiroga et al. (2004)
%     cfg.pISI = 1.5; % min ISI interval in ms  %peakseak Peter O'Connor
%     tdt_extractspikes_AP(cfg);
    
    %% ESIload the files
    pathLfp = dir(fullfile(outputDirSes, '*_xWav.lfp'));
%     pathMua = dir(fullfile(outputDirSes, '*_xWav.mua'));
%     pathMuax = dir(fullfile(outputDirSes, '*_xWav.muax'));
%     pathMuat = dir(fullfile(outputDirSes, '*_xWav.muat'));
    ESIload(fullfile(outputDirSes, pathLfp.name), '-mat');
    lfp.data = data; clear data;
%     ESIload(fullfile(outputDirSes, pathMuax.name), '-mat');
%     muax.data = data; clear data;
%     ESIload(fullfile(outputDirSes, pathMuat.name), '-mat');
%     muat.data = data; clear data;
%     ESIload(fullfile(outputDirSes, pathMua.name), '-mat');
%     mua.data = data; clear data;
    % ESIload TDT strobes
    sesName = sessionList{ses};
    savename = sprintf('%s/Experiments_%s', outputDirSes, sesName);
    [trialEnd, trialStart, trialStimOn, Cond, taccept] = getSnipsForSession(sesName, sesType);
    
    % Chop it - LFP
    data.fsample = lfp.data.fsample;
    data.cfg = lfp.data.cfg;
    data.label = lfp.data.label;
    data.trialinfo = [[1:length(Cond)]' taccept Cond trialStart trialEnd]; % add in trialno
    time = lfp.data.time{1};
    data.sampleinfo = [];
    for ii = 1:length(Cond)
        trialChosen = trialStart(ii) <= time & time <= trialEnd(ii);
        data.trial{ii} = lfp.data.trial{1}(:, trialChosen);
        data.time{ii} = time(trialChosen)-trialStimOn(ii);
        sinfo = round(data.fsample*[trialStart(ii) trialEnd(ii)]);
        sinfo(2) = sinfo(2) + length(data.time{ii}) - (sinfo(2)-sinfo(1))-1;
        data.sampleinfo = [data.sampleinfo; sinfo];
    end
    data.cfg.trl = data.sampleinfo;
    trialtime = [cellfun(@min, data.time)' cellfun(@max, data.time)'];
% %     save([savename '_chopped.lfp'], 'data', '-v7.3');
    clear data
    
    % filtered LFP
    cfg = [];
    cfg.bsfilter = 'yes';
    cfg.bsfreq   = [49.9 50.1;99.7 100.3;149.5 150.5;];
    cfg.bsinstabilityfix = 'reduce';
    dataFilt = ft_preprocessing(cfg, lfp.data);
    
    % check power spec
    cfg = [];
    cfg.length = 200;
    dataSnippets = ft_redefinetrial(cfg, dataFilt);
    
    cfg = [];
    cfg.foilim = [0 200];
    cfg.method = 'mtmfft';
    cfg.taper = 'rectwin';
    freq = ft_freqanalysis(cfg, dataSnippets);
    
    % plot the filtered spectrum
    h1 = figure(); set(h1, 'visible', 'off');
    for ii=1:32
        subplot(ceil(sqrt(32)), ceil(sqrt(32)), ii);
        loglog(freq.freq, abs(freq.powspctrm(ii, :)), 'r'); hold on;
        xlim(([1 max(freq.freq)]));
        ylim(([1e-16 max(abs(freq.powspctrm(ii,:)))]));
        grid on
        set(gca,'FontSize',5);
    end
    fname = [savename '_chopped_filtered.lfp.png'];
    print(h1, fname, '-dpng', '-r300');
    
    % Chop it - FilteredLFP
    data.fsample = lfp.data.fsample;
    data.cfg = lfp.data.cfg;
    data.label = lfp.data.label;
    data.trialinfo = [[1:length(Cond)]' taccept Cond trialStart trialEnd]; % add in trialno
    time = dataFilt.time{1};
    data.sampleinfo = [];
    for ii = 1:length(Cond)
        trialChosen = trialStart(ii) <= time & time <= trialEnd(ii);
        data.trial{ii} = dataFilt.trial{1}(:, trialChosen);
        data.time{ii} = time(trialChosen)-trialStimOn(ii);
        sinfo = round(data.fsample*[trialStart(ii) trialEnd(ii)]);
        sinfo(2) = sinfo(2) + length(data.time{ii}) - (sinfo(2)-sinfo(1))-1;
        data.sampleinfo = [data.sampleinfo; sinfo];
    end
    data.cfg.trl = data.sampleinfo;
    save([savename '_chopped_filtered.lfp'], 'data', '-v7.3');
    
% %     % Chop it - MUAX
% %     data.fsample = muax.data.fsample;
% %     data.cfg = muax.data.cfg;
% %     data.label = muax.data.label;
% %     data.trialinfo = [[1:length(Cond)]' taccept Cond trialStart trialEnd]; % add in trialno
% %     time = muax.data.time{1};
% %     data.sampleinfo = [];
% %     for ii = 1:length(Cond)
% %         trialChosen = trialStart(ii) <= time & time <= trialEnd(ii);
% %         data.trial{ii} = muax.data.trial{1}(:, trialChosen);
% %         data.time{ii} = time(trialChosen)-trialStimOn(ii);
% %         sinfo = round(data.fsample*[trialStart(ii) trialEnd(ii)]);
% %         sinfo(2) = sinfo(2) + length(data.time{ii}) - (sinfo(2)-sinfo(1))-1;
% %         data.sampleinfo = [data.sampleinfo; sinfo];
% %     end
% %     data.cfg.trl = data.sampleinfo;
% %     save([savename '_chopped.muax'], 'data', '-v7.3');
% %     clear data
% %     
% %     % Chop it - MUA
% %     data.fsample = mua.data.fsample;
% %     data.cfg = mua.data.cfg;
% %     data.label = mua.data.label;
% %     data.trialinfo = [[1:length(Cond)]' taccept Cond trialStart trialEnd]; % add in trialno
% %     time = mua.data.time{1};
% %     data.sampleinfo = [];
% %     for ii = 1:length(Cond)
% %         trialChosen = trialStart(ii) <= time & time <= trialEnd(ii);
% %         data.trial{ii} = mua.data.trial{1}(:, trialChosen);
% %         data.time{ii} = time(trialChosen)-trialStimOn(ii);
% %         sinfo = round(data.fsample*[trialStart(ii) trialEnd(ii)]);
% %         sinfo(2) = sinfo(2) + length(data.time{ii}) - (sinfo(2)-sinfo(1))-1;
% %         data.sampleinfo = [data.sampleinfo; sinfo];
% %     end
% %     data.cfg.trl = data.sampleinfo;
% %     save([savename '_chopped.mua'], 'data', '-v7.3');
% %     clear data
% %     
% %     % Chop it - MUAThreshold
% %     data.fsample = muat.data.fsample;
% %     data.cfg = muat.data.cfg;
% %     data.label = muat.data.label;
% %     data.trialinfo = [[1:length(Cond)]' taccept Cond trialStart trialEnd]; % add in trialno
% %     
% %     % Make Spikes
% %     spike.hdr = muat.data.hdr;
% %     spike.trialinfo = [[1:length(Cond)]' taccept Cond trialStart trialEnd]; % add in trialno
% %     nrTrials = length(Cond);
% %     data.sampleinfo = [];
% %     for ch = 1:length(muat.data.label)
% %         spikeTimes = muat.data.trial{ch};
% %         spikeStamps = round(data.fsample*muat.data.trial{ch});
% %         for ii = 1:nrTrials
% %             if ch == 1
% %                 data.sampleinfo = [data.sampleinfo; round(data.fsample*[trialStart(ii) trialEnd(ii)])];
% %             end
% %             data.trial{ch, ii} = spikeTimes(and(spikeTimes>trialStart(ii), spikeTimes<=trialEnd(ii)))- trialStimOn(ii);
% %             data.timestamp{ch, ii} = spikeStamps(and(spikeTimes>trialStart(ii), spikeTimes<=trialEnd(ii)))- trialStimOn(ii);
% %         end
% %         data.cfg.trl = data.sampleinfo;
% %         spike.time{ch} = ([data.trial{ch, :}]');
% %         spike.timestamp{ch} = ([data.timestamp{ch, :}]');
% %         spike.label{ch} = data.label{ch};
% %         spike.trialtime = trialtime;
% %         spike.cond = Cond;
% %         spike.sampleinfo = data.sampleinfo;
% %         spike.taccept = taccept;
% %         onesVec = cellfun(@(x) (ones(length(x), 1)'), data.trial(ch, :), 'UniformOutput', false);
% %         trialVec = cellfun(@(x, y) (x*y), onesVec, num2cell([1:nrTrials]), 'UniformOutput', false);
% %         spike.trial{ch} = [trialVec{:}]';
% %     end
% %     save([savename '_chopped.muat'],'data', '-v7.3')
% %     clear mua
% %     
% %     save([savename '_chopped.spike'],'spike', '-v7.3')
% %     clear spike
    
    %     %% Get the eye data
    %     if str2num(sesName(6:8)) > 458 | strcmp(sesName(1:4), 'ares')
    %         % Convert the EDF file to ASC
    %         pathEDF = dir(fullfile(outputDirSes, '*.edf'));
    %         system(sprintf('linuxinstall/bin/edf2asc -t -c %s', fullfile(outputDirSes, pathEDF.name)));
    %         %  -sp
    %         % Get the ASC
    %         eyeList = dir(fullfile(outputDirSes, '*Eye.asc'));
    %         cfg = [];
    %         cfg.dataset = fullfile(outputDirSes, eyeList.name);
    %         cfg.trialOutcome = taccept;
    %         cfg.savepath = savename;
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
    %         data.trialinfo = [[1:length(Cond)]' taccept Cond trialStart trialEnd]; % add in trialno
    %         %     time  = 0:1/data.fsample:(1/data.fsample)*(size(Eyes, 2)-1);
    %         time = dataEye.time{1};
    %         for ii = 1:length(Cond)
    %             trialChosen = trialStart(ii) <= time & time <= trialEnd(ii);
    %             data.trial{ii} = dataEye.trial{1}(:, trialChosen);
    %             data.time{ii} = time(trialChosen)-trialStart(ii)-preStim;
    %         end
    %         ESIsave([savename 'StimOn.eyes'],'data');
    %         clear data
    %     end
    %     catch
    %         failedSessions = [failedSessions; {sesName}]
    %     end
end
% failedSessions
