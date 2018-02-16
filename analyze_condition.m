function analyze_condition(allCfg, lfpClean, muaClean, spikeClean, trialsClean, trialsChosen, varargin)
% Primer
if nargin==8
    cnd = varargin{1};
    cnf = varargin{2};
    fbase = sprintf('cond_%02d_%.2f_%s', cnd, cnf, allCfg.tag);
else
    cnd = varargin{1};
    fbase = sprintf('cond%02d_%s', cnd, allCfg.tag);
end
nChan = length(lfpClean.label);

%% Create LFP MUA Trials for full trial
cfg = [];
cfg.trials = trialsChosen;
cfg.toilim = allCfg.timeAll;
lfpSel = ft_redefinetrial(cfg, lfpClean);

cfg = [];
cfg.trials = trialsChosen;
cfg.toilim = allCfg.timeAll;
muaSel = ft_redefinetrial(cfg, muaClean);

%% run stuff
if allCfg.runTimelockLFP; runTimelockLFP(allCfg, fbase, lfpSel); end
if allCfg.runTimelockMUAX; runTimelockMUAX(allCfg, fbase, muaSel); end
if allCfg.runPSTH; runPSTH(allCfg, fbase, trialsClean, trialsChosen, spikeClean); end

%% Create LFP MUA Trials for sustained trial
cfg = [];
cfg.trials = trialsChosen;
cfg.toilim = allCfg.timeSustained;
lfpSel = ft_redefinetrial(cfg, lfpClean);

% Get event
cfg = [];
cfg.trials = trialsChosen;
cfg.toilim = allCfg.timeSustained;
muaSel = ft_redefinetrial(cfg, muaClean);

%% Computes the multitaper fft
if allCfg.runTFR; runTFR(allCfg, fbase, trialsChosen, lfpSel); end
if allCfg.runConnectivity; runConnectivity(allCfg, fbase, trialsChosen, lfpSel, muaSel); end
% if allCfg.runSFC; runSFC(spikeClean, lfpSel); end
% if allCfg.runSTA; runSTA(spikeClean, lfpSel); end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%           Subfunctions           %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Timelock LFP
function runTimelockLFP(allCfg, fbase, lfpSel)
cfg = [];
cfg.keeptrials = 'no';
cfg.removemean = 'yes';
cfg.vartrllength = 2;
timelockLFP = ft_timelockanalysis(cfg, lfpSel);

cfg = [];
cfg.baseline = allCfg.timeBaseline;
timelockLFP = ft_timelockbaseline(cfg, timelockLFP);

% save
fname = sprintf('%stimelockLFP.mat', fbase);
ESIsave(fullfile(allCfg.outputfile, fname), 'timelockLFP');
clear timelockLFP
return

%% Timelock MUAX
function runTimelockMUAX(allCfg, fbase, muaSel)
cfg = [];
cfg.keeptrials = 'no';
cfg.removemean = 'yes';
cfg.vartrllength = 2;
timelockMUAX = ft_timelockanalysis(cfg, muaSel);

cfg = [];
cfg.baseline = allCfg.timeBaseline;
timelockMUAX = ft_timelockbaseline(cfg, timelockMUAX);

cfg = [];
cfg.timwin = [-0.005 0.005];
cfg.spikechannel = muaSel.label;
timelockMUAX = ft_spikedensity(cfg, timelockMUAX);

%save
fname = sprintf('%stimelockMUAX.mat', fbase);
ESIsave(fullfile(allCfg.outputfile, fname), 'timelockMUAX');
clear timelockMUAX
return

%% PSTH
function runPSTH(allCfg, fbase, trialsClean, trialsChosen, spikeClean)
spikeClean.trialtime(:, 1) = allCfg.timeAll(1);
spikeClean.trialtime(:, 2) = allCfg.timeAll(2);

cfg = [];
cfg.binsize = 0.02;
cfg.trials = [trialsClean(trialsChosen)];
cfg.keeptrials = 'no';
trialPSTH = ft_spike_psth(cfg, spikeClean);

cfg = [];
cfg.baseline = [-0.5 0];
cfg.keeptrials = 'yes';
trialPSTH = ft_timelockbaseline(cfg, trialPSTH);
% 
% cfg = [];
% cfg.timwin = [-0.005 0.005];
% cfg.spikechannel = spikeClean.label;
% trialPSTH = ft_spikedensity(cfg, trialPSTH);

% save these
fname = sprintf('%strialPSTH.mat', fbase);
ESIsave(fullfile(allCfg.outputfile, fname), 'trialPSTH');
clear trialPSTH

if allCfg.runErrorBars
    cfg = [];
    cfg.binsize = 0.02;
    cfg.trials = trialsClean(trialsChosen);
    cfg.keeptrials = 'no';
    
    trialPSTHVar = ft_spike_psth(cfg, spikeClean);
    fname = sprintf('%strialPSTHVar.mat', fbase);
    ESIsave(fullfile(allCfg.outputfile, fname), 'trialPSTHVar');
end
return

%% TFR
function runTFR(allCfg, fbase, trialsChosen, lfpSel)
if strcmp(allCfg.name, 'HermesAA'); lfpSel = rmfield(lfpSel, 'label_tdt'); end; % unfortunately
cfg = [];
cfg.output     = 'pow';
cfg.channel    = 'all';
cfg.method     = 'mtmfft';
if strcmp(allCfg.type, 'NatImSEQ')
    cfg.foi        = 10:2.5:120;
else
    cfg.foi        = 10:2:120;
end
cfg.channel = lfpSel.label;
cfg.taper      = 'dpss';
cfg.tapsmofrq = 7*ones(1,length(cfg.foi));

lfpPower = ft_freqanalysis(cfg, lfpSel);
if ~allCfg.flagSecond
    fname = sprintf('%slfpPower.mat', fbase);
    ESIsave(fullfile(allCfg.outputfile, fname), 'lfpPower');
    clear lfpPower
else
    lfpPower_Second = lfpPower;
    fname = sprintf('%slfpPower_Second.mat', fbase);
    ESIsave(fullfile(allCfg.outputfile, fname), 'lfpPower_Second');
    clear lfpPower_Second
end

if allCfg.runErrorBars
    % Computes the multitaper fft EVENT
    cfg = [];
    cfg.output     = 'pow';
    cfg.channel    = 'all';
    cfg.method     = 'mtmfft';
    cfg.foi        = 10:2:120;
    cfg.channel = lfpSel.label;
    cfg.taper      = 'dpss';
    cfg.tapsmofrq = 7*ones(1,length(cfg.foi));
    
    for tr=1:length(trialsChosen)
        cfg.trials = ~ismember(trialsChosen, trialsChosen(tr));
        lfpPowerVar(tr) = ft_freqanalysis(cfg, lfpSel);
    end
    if ~allCfg.flagSecond
        fname = sprintf('%slfpPowerVar.mat', fbase);
        ESIsave(fullfile(allCfg.outputfile, fname), 'lfpPowerVar');
        clear lfpPowerVar
    else
        lfpPowerVar_Second = lfpPowerVar;
        fname = sprintf('%slfpPowerVar_Second.mat', fbase);
        ESIsave(fullfile(allCfg.outputfile, fname), 'lfpPowerVar_Second');
        clear lfpPowerVar_Second
    end
    clear TFRmultEvent
end
return

%% Connectivity
function runConnectivity(allCfg, fbase, trialsChosen, lfpSelAll, muaSelAll)
muaSelAll.label = strcat('muax-', muaSelAll.label);
if strcmp(allCfg.name, 'Hermes') && isfield(muaSelAll, 'label_tdt')
    label_tdt = muaSelAll.label_tdt;
    muaSelAll = rmfield(muaSelAll, 'label_tdt');
    lfpSelAll = rmfield(lfpSelAll, 'label_tdt');
end % unfortunately
%         for ch=1:length(muaSelAll.label)
%             cfg = [];
%             cfg.channel = muaSelAll.label(ch);
%             muaSelCh = ft_selectdata(cfg, muaSelAll);
adata = ft_appenddata([], lfpSelAll, muaSelAll);

for ch=1:length(muaSelAll.label)
    % get spectrum
    cfg            = [];
    cfg.output     = 'powandcsd';
    cfg.method     = 'mtmfft';
    cfg.foi        = 10:2:120;
    cfg.taper      = 'dpss';
    cfg.tapsmofrq  = 7*ones(1,length(cfg.foi));
    cfg.keeptrials = 'yes';
    cfg.channel    = 'all';
    cfg.channelcmb = [repmat(muaSelAll.label(ch), length(lfpSelAll.label), 1) lfpSelAll.label];
    freq           = ft_freqanalysis(cfg, adata);
    
    % get stats
    cfg            = [];
    cfg.channelcmb = [repmat(muaSelAll.label(ch), length(lfpSelAll.label), 1) lfpSelAll.label];
    cfg.method     = 'coh';
    muaxCoherence(ch)  = ft_connectivityanalysis(cfg, freq);
end
if ~allCfg.flagSecond
    fname = sprintf('%smuaxCoherence.mat', fbase);
    ESIsave(fullfile(allCfg.outputfile, fname), 'muaxCoherence');
    clear muaxCoherence
else
    muaxCoherence_Second = muaxCoherence;
    fname = sprintf('%smuaxCoherence_Second.mat', fbase);
    ESIsave(fullfile(allCfg.outputfile, fname), 'muaxCoherence_Second');
    clear muaxCoherence_Second
end
return
