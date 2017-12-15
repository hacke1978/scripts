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
savename = allCfg.outputfile

% Load the data
if strcmp(allCfg.name, 'hermes')
    pathLFP = dir(fullfile(filename, '*stimOn.lfp'));
    pathMUAX = dir(fullfile(filename, '*stimOn.muax'));
    pathSpike = dir(fullfile(filename, sprintf('*_%s_chopped.spike', allCfg.tag)));
else
    pathLFP = dir(fullfile(filename, '*_chopped.lfp'));
    pathMUAX = dir(fullfile(filename, '*_chopped.muax'));
    pathSpike = dir(fullfile(filename, sprintf('*_%s_chopped.spike', allCfg.tag)));
end
load(fullfile(filename, pathLFP.name), '-mat');
lfp.data = data; clear data;
load(fullfile(filename, pathMUAX.name), '-mat');
muax.data = data; clear data;
load(fullfile(filename, pathSpike.name), '-mat');

% get conditions / channels
if strcmp(allCfg.name, 'Hermes')
    if strcmp(allCfg.type, 'grating-ori')
        taccept = (muax.data.trialinfo(:, 4)==1);
        CondSF = muax.data.trialinfo(:, 3);
    elseif strcmp(allCfg.type, 'rfmapping-bar')
        taccept = (ones([size(muax.data.trialinfo, 1) 1]));
    else
        taccept = (muax.data.trialinfo(:, 8)==1);
    end
    v1accept = strncmp(muax.data.label, 'V1', 2);
    caccept = v1accept;
    Cond = muax.data.trialinfo(:, 2);
else
    taccept = (muax.data.trialinfo(:, 2)==0);
    if strcmp(allCfg.name, 'Ares')
        caccept = ones(size(muax.data.label));
        caccept(cellfun(@(x) (str2num(x(8:end))>32), muax.data.label)) = false;
        caccept = logical(caccept);
    else
        caccept = logical(ones(size(muax.data.label)));
    end
    Cond = muax.data.trialinfo(:, 3);
end
%% Preprocessing
% Remove 50HzLine Noise from RAW LFP Data
if allCfg.filterLineNoise
    if isempty(dir(fullfile(filename, '*_choppedbsline.lfp')))
        
        % get full trial
        pathLFP = dir(fullfile(filename, '*xWav.lfp'));
        lfpfull = load(fullfile(filename, pathLFP.name), '-mat');
        
        % remove line noise
        cfg = [];
        cfg.bsfilter = 'yes';
        cfg.bsfreq   = [49.9 50.1];
        cfg.bsinstabilityfix = 'reduce';
        lfpfull.data = ft_preprocessing(cfg, lfpfull.data);
        
        %Chop into trials
        tok = strsplit(filename, '/');
        [trialEnd, trialStart, trialStimOn, Cond, taccept] = getSnipsForSession(tok{end}, allCfg.type);
        taccept = taccept == 0;
        clear lfp
        data  = chop_it(lfpfull, trialEnd, trialStart, Cond, taccept, trialStimOn);
        
        saveName = fullfile(filename, pathLFP.name(1:end-4));
        save([saveName '_choppedbsline.lfp'],'data', '-v7.3');
        clear data lfpfull
    end
    % load filtered lfp
    pathLFP = dir(fullfile(filename, '*_choppedbsline.lfp'));
    lfp = load(fullfile(filename, pathLFP.name), '-mat');
end

if allCfg.getBadChannels
    % get full trial
    pathLFP = dir(fullfile(filename, '*xWav.lfp'));
    lfpfull = load(fullfile(filename, pathLFP.name), '-mat');
    
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

if strcmp(allCfg.type, 'grating-ori')
    CondSF = CondSF(taccept);
    condLibSF = unique(CondSF);
    for sf=1:length(condLibSF)
        for cnd=1:length(condLib)
            cnd
            trialsChosen = find(Cond == condLib(cnd) & CondSF == condLibSF(sf));
            tok = strsplit(filename, '/');
            if strcmp(tok{end}, 'hermes_20170428_fixation-naturalim_15')
                trialsChosen = trialsChosen(1:18);
            end
            analyze_condition(allCfg, lfpClean, muaClean, spikeClean, trialsClean, trialsChosen, condLib(cnd), condLibSF(sf));
        end
    end
elseif strcmp(allCfg.type, 'NatImSEQ')
    for cnd=1:length(condLib)
        cnd
        trialsChosen = find(Cond==condLib(cnd));
        %first stim
        allCfg.flagSecond = false;
        allCfg.timeSustained = allCfg.timeSustained_First;
        analyze_condition(allCfg, lfpClean, muaClean, spikeClean, trialsClean, trialsChosen, condLib(cnd))
        
        % second stim
        allCfg.flagSecond = true;
        allCfg.timeSustained = allCfg.timeSustained_Second;
        analyze_condition(allCfg, lfpClean, muaClean, spikeClean, trialsClean, trialsChosen, condLib(cnd))
    end
else
    for cnd=1:length(condLib)
        trialsChosen = find(Cond==condLib(cnd))';
        length(trialsChosen)
        tok = strsplit(filename, '/');
        if strcmp(tok{end}, 'hermes_20170428_fixation-naturalim_15')
            trialsChosen = trialsChosen(1:18);
        end
        analyze_condition(allCfg, lfpClean, muaClean, spikeClean, trialsClean, trialsChosen, condLib(cnd))
    end
end

% Save Baseline TFR
if (allCfg.runBaseline && allCfg.runTFR); runTFR_Baseline(); end
if (allCfg.runSFC && allCfg.Baseline); runSFC_Baseline(); end
if (allCfg.STA && allCfg.Baseline); runSTA_Baseline(); end
end
%% Subfunctions

%% TFR Baseline
function runTFR_Baseline()
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
function runSFC_Baseline()
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
cfg.method = 'mtmconvol';
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
function runSTA_Baseline()
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
ESIsave(fullfile(outputDirSes, 'data', sprintf('im_%02d_staBSEvent.mat', imId)), 'staBSEvent');
end