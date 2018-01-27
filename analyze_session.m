function analyze_session(allCfg)
addpath('/mnt/hpx/slurm/uranc/fieldtrip/');
ft_defaults

if ~isfield(allCfg, 'filterLineNoise');    allCfg.filterLineNoise = false; end
if ~isfield(allCfg, 'getBadChannels');     allCfg.getBadChannels = false;  end
if ~isfield(allCfg, 'runPSTH');            allCfg.timelockLfp = false;     end
if ~isfield(allCfg, 'runTimelockLFP');     allCfg.timelockLfp = false;     end
if ~isfield(allCfg, 'runTimelockMUAX');    allCfg.timelockMuax = false;    end
if ~isfield(allCfg, 'runTFR');             allCfg.tfr = false;             end
if ~isfield(allCfg, 'save');               allCfg.save = false;            end

% Load neural data - LFP, MUA, Waveform

filename = allCfg.inputfile;
if ~iscell(filename); filename = {filename}; end;

% merge sessions if any
if iscell(allCfg.outputfile)
    if strcmp(allCfg.name, 'Hermes')
        tok = cellfun(@(x) strsplit(x, '_'), allCfg.outputfile, 'UniformOutput', false);
        tok = vertcat(tok{:});
        mergeName  = allCfg.outputfile{1};
        for ii=2:length(allCfg.outputfile)
            mergeName = [mergeName , sprintf('_%s', tok{ii, end})];
        end
    else
        tok = cellfun(@(x) strsplit(x, '/'), allCfg.outputfile, 'UniformOutput', false);
        tok = vertcat(tok{:});
        mergeName  = allCfg.outputfile{1};
        for ii=2:length(allCfg.outputfile)
            mergeName = [mergeName , sprintf('_%s', tok{ii, end}(end-2:end))];
        end
    end
    allCfg.outputfile = mergeName;
end

% saveDir
if ~exist(allCfg.outputfile, 'dir')
    mkdir(allCfg.outputfile);
end

% Load the data - handles session merging
for ii=1:length(filename);
    % Remove 50HzLine Noise from RAW LFP Data
    %     if allCfg.filterLineNoise
    %         if isempty(dir(fullfile(filename{ii}, '*_choppedbsline.lfp')))
    %
    %             % get full trial
    %             pathLFP = dir(fullfile(filename{ii}, '*xWav.lfp'));
    %             lfpfull = load(fullfile(filename{ii}, pathLFP.name), '-mat');
    %
    %             % remove line noise
    %             %% USE Same
    %             cfg = [];
    %             cfg.bsfilter = 'yes';
    %             cfg.bsfreq   = [49.9 50.1];
    %             cfg.bsinstabilityfix = 'reduce';
    %             lfpfull.data = ft_preprocessing(cfg, lfpfull.data);
    %
    %             %Chop into trials
    %             tok = strsplit(filename, '/');
    %             [trialEnd, trialStart, trialStimOn, Cond, taccept] = getSnipsForSession(tok{end}, allCfg.type);
    %             taccept = taccept == 0;
    %             clear lfp
    %             data  = chop_it(lfpfull, trialEnd, trialStart, Cond, taccept, trialStimOn);
    %
    %             saveName = fullfile(filename{ii}, pathLFP.name(1:end-4));
    %             save([saveName '_choppedbsline.lfp'],'data', '-v7.3');
    %             clear data lfpfull
    %         end
    %     end
    
    % get filenamaes
    if allCfg.filterLineNoise
        pathLFP = dir(fullfile(filename{ii}, '*_chopped_filtered.lfp'));
    else
        pathLFP = dir(fullfile(filename{ii}, '*_chopped.lfp'));
    end
    pathMUAX = dir(fullfile(filename{ii}, '*_chopped.muax'));
    pathSpike = dir(fullfile(filename{ii}, sprintf('*%s_chopped.spike', allCfg.tag)));
    
    if ii==1
        load(fullfile(filename{ii}, pathLFP.name), '-mat');
        lfp.data = data; clear data;
        load(fullfile(filename{ii}, pathMUAX.name), '-mat');
        muax.data = data; clear data;
        load(fullfile(filename{ii}, pathSpike.name), '-mat');
    else
        load(fullfile(filename{ii}, pathLFP.name), '-mat');
        lfp.data = ft_appenddata([], lfp.data, data);
        load(fullfile(filename{ii}, pathMUAX.name), '-mat');
        muax.data = ft_appenddata([], muax.data, data);
        
        % spikes are a little more complicated
        spikeAdd = load(fullfile(filename{ii}, pathSpike.name), '-mat');
        spike.trialinfo = [spike.trialinfo; spikeAdd.spike.trialinfo];
        spike.sampleinfo = [spike.sampleinfo; spikeAdd.spike.sampleinfo];
        spike.trialtime = [spike.trialtime; spikeAdd.spike.trialtime];
        spike.cond = [spike.cond; spikeAdd.spike.cond];
        spike.taccept = [spike.taccept; spikeAdd.spike.taccept];
        nTrials = max(spike.trial{1});
        for ch = 1:length(spike.label)
            spike.time{ch} = [spike.time{ch}; spikeAdd.spike.time{ch}];
            spike.timestamp{ch} = [spike.timestamp{ch}; spikeAdd.spike.timestamp{ch}];
            spike.trial{ch} = [spike.trial{ch}; spikeAdd.spike.trial{ch}+nTrials];
        end
    end
end

% if length(filename)>1
%     % get merge dirs
%     tok = strfind(filename{1}, '/');
%     mname = strsplit(mergeName, '/');
%     mkdir(fullfile(filename{1}(1:tok(end)), mname{end}));
%     mergeDir = fullfile(filename{1}(1:tok(end)), mname{end}, mname{end});
%
%     % save lfp
%     data = lfp.data;
%     save([mergeDir '_chopped.lfp'],'data', '-v7.3');
%     clear data;
%     % save muax
%     data = muax.data;
%     save([mergeDir '_chopped.muax'],'data', '-v7.3');
%     clear data;
%     % save spike
%     save([mergeDir '_chopped.spike'],'spike', '-v7.3');
% end

% get conditions / channels
if strcmp(allCfg.name, 'Hermes')
    if strcmp(allCfg.type, 'grating-ori')
        taccept = (muax.data.trialinfo(:, 4)==1);
        CondSF = muax.data.trialinfo(:, 3);
    elseif strcmp(allCfg.type, 'rfmapping-bar')
        taccept = (ones([size(muax.data.trialinfo, 1) 1]));
    else
        %             taccept = (muax.data.trialinfo(:, 8)==1);
        taccept = (muax.data.trialinfo(:, 2)==1);
    end
    v1accept = strncmp(muax.data.label, 'V1', 2);
    caccept = v1accept;
    %         Cond = muax.data.trialinfo(:, 2);
    Cond = muax.data.trialinfo(:, 3);
else
    taccept = (muax.data.trialinfo(:, 2)==0);
    if strcmp(allCfg.name, 'Ares')
        caccept = ones(size(muax.data.label));
        caccept(cellfun(@(x) (str2num(x(6:end))>32), muax.data.label)) = false;
        caccept = logical(caccept);
    else
        caccept = logical(ones(size(muax.data.label)));
    end
    Cond = muax.data.trialinfo(:, 3);
end
%% Check for bad channels

if allCfg.getBadChannels
    % get full trial
    pathLFP = dir(fullfile(filename{ii}, '*xWav.lfp'));
    lfpfull = load(fullfile(filename{ii}, pathLFP.name), '-mat');
    
    % check for bad channels
    cfg = [];
    cfg.foilim = [0 200];
    cfg.method = 'mtmfft';
    cfg.taper = 'rectwin';
    cfg.pad = 'nextpow2';
    freq = ft_freqanalysis(cfg, lfpfull.data);
    
    % check for bad channels
    fsel = nearest(freq.freq, 50);
    z = zscore(freq.powspctrm(:,fsel));
    caccept = caccept & abs(z)<3;
    sprintf('Number of bad channels: %d', sum(~caccept));
end

%% Get Clean Trials

%selecdct
cfg = [];
cfg.trials = taccept';
cfg.nanmean = 'yes';
cfg.channel = find(caccept);
muaClean = ft_selectdata(cfg, muax.data);

cfg = [];
cfg.trials = taccept';
cfg.channel = find(caccept);
lfpClean = ft_selectdata(cfg, lfp.data);


cfg = [];
cfg.trials = find(taccept)';
cfg.spikechannel = find(caccept)';
spikeClean  = ft_spike_select(cfg, spike);

%     cfg = [];
%     cfg.trials = taccept;
%     cfg.channel = find(caccept);
%     %     cfg.toilim = [-0.5 1.2];
%     eyeClean = ft_selectdata(cfg, eye.data);
%
%     cfg = [];
%     cfg.trials = taccept;
%     cfg.channel = find(caccept);
%     %     cfg.toilim = [-0.5 1.2];
%     pdClean = ft_selectdata(cfg, pd.data);

%% Get all conditions
Cond = Cond(taccept);
trialsClean = find(taccept);
condLib = unique(Cond);

% % if strcmp(allCfg.type, 'grating-ori')
% %     CondSF = CondSF(taccept);
% %     condLibSF = unique(CondSF);
% %     for sf=1:length(condLibSF)
% %         for cnd=1:length(condLib)
% %             cnd
% %             trialsChosen = find(Cond == condLib(cnd) & CondSF == condLibSF(sf));
% %             %                 tok = strsplit(filename, '/');
% %             %                 if strcmp(tok{end}, 'hermes_20170428_fixation-naturalim_15')
% %             %                     trialsChosen = trialsChosen(1:18);
% %             %                 end
% %             analyze_condition(allCfg, lfpClean, muaClean, spikeClean, trialsClean, trialsChosen, condLib(cnd), condLibSF(sf));
% %         end
% %     end
% % elseif strcmp(allCfg.type, 'NatImSEQ')
% %     %     maskFirst = muax.data.trialinfo(taccept, end-5);
% %     %     maskSize = muax.data.trialinfo(taccept, end-4);
% %     %     maskPos = muax.data.trialinfo(taccept, end-2);
% %     %     maskSize(maskSize==6 & maskPos==2) = 6.01;
% %     %     maskLib = unique(maskSize);
% %     %     seqLib = unique(maskFirst);
% %     for cnd=1:length(condLib)
% %         %         for ms=1:length(maskLib)
% %         %             for sq=1:length(seqLib)
% %         %                 trialsChosen = find(Cond==condLib(cnd) & maskSize == maskLib(ms) & maskFirst == seqLib(sq));
% %         trialsChosen = find(Cond==condLib(cnd));
% %         %first stim
% %         allCfg.flagSecond = false;
% %         allCfg.timeSustained = allCfg.timeSustained_First;
% %         analyze_condition(allCfg, lfpClean, muaClean, spikeClean, trialsClean, trialsChosen, condLib(cnd));
% %         %                     seqLib(sq)*(length(condLib)*length(maskLib))+length(maskLib)*(cnd-1)+ms)
% %         
% %         % second stim
% %         allCfg.flagSecond = true;
% %         allCfg.timeSustained = allCfg.timeSustained_Second;
% %         analyze_condition(allCfg, lfpClean, muaClean, spikeClean, trialsClean, trialsChosen, condLib(cnd));
% %         %                     seqLib(sq)*(length(condLib)*length(maskLib))+length(maskLib)*(cnd-1)+ms)
% %         %             end
% %         %         end
% %     end
% % else
% %     for cnd=1:length(condLib)
% %         trialsChosen = find(Cond==condLib(cnd))';
% %         length(trialsChosen)
% %         %             tok = strsplit(filename, '/');
% %         %             if strcmp(tok{end}, 'hermes_20170428_fixation-naturalim_15')
% %         %                 trialsChosen = trialsChosen(1:18);
% %         %             end
% %         analyze_condition(allCfg, lfpClean, muaClean, spikeClean, trialsClean, trialsChosen, condLib(cnd))
% %     end
% % end

% Save Baseline TFR
% if (allCfg.runTFR && allCfg.runBaseline); runTFR_Baseline(allCfg, lfpClean); end
if (allCfg.runSFC && allCfg.runBaseline); runSFC_Baseline(allCfg, lfpClean, spikeClean); end
if (allCfg.runSTA && allCfg.Baseline); runSTA_Baseline(allCfg, lfpClean, spikeClean); end
end

%%% Subfunctions
%% TFR Baseline
function runTFR_Baseline(allCfg, lfpClean)
savename = allCfg.outputfile;

cfg = [];
cfg.trials = 'all';
cfg.toilim = allCfg.timeBaseline
lfpBaseline = ft_redefinetrial(cfg, lfpClean);

cfg = [];
cfg.output     = 'pow';
cfg.channel = lfpBaseline.label;
cfg.method     = 'mtmfft';
if strcmp(allCfg.type, 'NatImSEQ')
    cfg.foi        = 10:2.5:120;
else
    cfg.foi        = 10:2:120;
end
cfg.taper      = 'dpss';
cfg.tapsmofrq = 7*ones(1,length(cfg.foi));
powerBaseline = ft_freqanalysis(cfg, lfpBaseline);
ESIsave(fullfile(savename, sprintf('cond_All_powerBaseline.mat')), 'powerBaseline');
clear powerBaseline
end

%% SFC Baseline
function runSFC_Baseline(allCfg, lfpClean, spikeClean)
savename = allCfg.outputfile;

cfg = [];
cfg.trials = 'all';
cfg.toilim = allCfg.timeBaseline
lfpBaseline = ft_redefinetrial(cfg, lfpClean);

% Select spike
nChan = length(spikeClean.trial);
cfg = [];
cfg.toilim = allCfg.timeBaseline;
spikeBaseline = ft_spike_select(cfg, spikeClean);

trl = unique([spikeBaseline.trial{:}]);
for ch = 1:nChan
    for tr = 1:length(trl)
        spikeBaseline.trial{ch}(spikeBaseline.trial{ch}==trl(tr)) = tr;
    end
end
spikeBaseline.label = strcat('label ',spikeBaseline.label);
spike.trialtime = repmat(allCfg.timeBaseline, length(trl), 1);
spikeBaseline.timestamp = spikeBaseline.time;

% % % % Baseline SFC
cfg = [];
% cfg.method = 'mtmconvol';
cfg.method = 'mtmfft';
cfg.foi    = 2:2:120;
cfg.t_ftimwin = 7./cfg.foi; % watch out that the first frequency is not too long duration, so have to set to a minimum value
twin = allCfg.timeSustained(2)-allCfg.timeSustained(1);
cfg.t_ftimwin(cfg.t_ftimwin>twin) = twin;
cfg.taper     = 'hanning';
[stsBase] = ft_spiketriggeredspectrum(cfg, lfpBaseline, spikeBaseline);

for k=1:nChan
    cfg = [];
    cfg.method = 'ppc1'; % use ppc1
    cfg.spikechannel = k;
    stSpecBaseline(k) = ft_spiketriggeredspectrum_stat(cfg, stsBase);
end

% Save baseline
ESIsave(fullfile(savename, 'cond_All_stSpecBaseline.mat'), 'stSpecBaseline');
clear stSpecBaseline
end

%% STA Baseline
function runSTA_Baseline(allCfg, lfpClean, spikeClean)
savename = allCfg.outputfile;

% Append spike
data_lfp = baselineDataLP;
data_lfp.hdr = data.hdr;
baselineSpike.timestamp = baselineSpike.time;
baselineSpike.trialtime = baselineSpike.basetime;
data_lfp.cfg.trl = baselineSpike.basetime;
data_all_base = ft_appendspike([], data_lfp, baselineSpike);

% Spike triggered average
clear staPre
cfg              = [];
cfg.keeptrials   = 'no';
cfg.latency      = [-0.5 0];
for k = 1:32
    if accept(k)
        cfg.spikechannel = baselineSpike.label{k}; % first unit
        cfg.channel      = data_lfp.label(k); % first four chans
        cfg.timwin = [-0.2 0.2]; % take 400 ms
        staBSEvent(k)        = ft_spiketriggeredaverage(cfg, data_all_base);
    end
end
% Save baseline STA
ESIsave(fullfile(allCfg.outputfile, 'data', sprintf('im_%02d_staBSEvent.mat', imId)), 'staBSEvent');
end