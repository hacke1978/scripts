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

savefile = allCfg.outputfile;
nChan = length(lfpClean.label);

% Create LFP MUA Trials for full trial
cfg = [];
cfg.trials = trialsChosen;
cfg.toilim = allCfg.timeAll;
lfpSel = ft_redefinetrial(cfg, lfpClean);

cfg = [];
cfg.trials = trialsChosen;
cfg.toilim = allCfg.timeAll;
muaSel = ft_redefinetrial(cfg, muaClean);

% run stuff
if allCfg.runTimelockLFP; runTimelockLFP(lfpSel); end
if allCfg.runTimelockMUAX; runTimelockMUAX(muaSel); end
if allCfg.runPSTH; runPSTH(spikeClean); end

% Create LFP MUA Trials for sustained trial
cfg = [];
cfg.trials = trialsChosen;
cfg.toilim = allCfg.timeSustained;
lfpSel = ft_redefinetrial(cfg, lfpClean);

% Get event
cfg = [];
cfg.trials = trialsChosen;
cfg.toilim = allCfg.timeSustained;
muaSel = ft_redefinetrial(cfg, muaClean);

% Computes the multitaper fft
if allCfg.runTFR; runTFR(lfpSel); end
if allCfg.runConnectivity; runConnectivity(lfpSel, muaSel); end
if allCfg.runSFC; runSFC(spikeClean, lfpSel); end
if allCfg.runSTA; runSTA(spikeClean, lfpSel); end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%           Subfunctions           %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Timelock LFP
    function runTimelockLFP(lfpSel)
        cfg = [];
        cfg.keeptrials = 'no';
        cfg.removemean = 'yes';
        cfg.vartrllength = 2;
        timelockLFP = ft_timelockanalysis(cfg, lfpSel);
        
        cfg = [];
        cfg.baseline = allCfg.timeBaseline;
        timelockLFP = ft_timelockbaseline(cfg, timelockLFP);
        
        % save
        fname = sprintf('%s_timelockLFP.mat', fbase);
        ESIsave(fullfile(savefile, fname), 'timelockLFP');
        clear timelockLFP
    end

%% Timelock MUAX
    function runTimelockMUAX(lfpSel)
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
        fname = sprintf('%s_timelockMUAX.mat', fbase);
        ESIsave(fullfile(savefile, fname), 'timelockMUAX');
        clear timelockMUAX
    end

%% PSTH
    function runPSTH(spikeClean)
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
        
        % save these
        fname = sprintf('%s_trialPSTH.mat', fbase);
        ESIsave(fullfile(savefile, fname), 'trialPSTH');
        clear trialPSTH
        
        if allCfg.runErrorBars
            cfg = [];
            cfg.binsize = 0.02;
            cfg.trials = trialsClean(trialsChosen);
            cfg.keeptrials = 'no';
            
            trialPSTHVar = ft_spike_psth(cfg, spikeClean);
            fname = sprintf('%s_trialPSTHVar.mat', fbase);
            ESIsave(fullfile(savefile, fname), 'trialPSTHVar');
        end
    end

%% TFR
    function runTFR(lfpSel)
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
        cfg.tapsmofrq = 7;
        
        lfpPower = ft_freqanalysis(cfg, lfpSel);
        fname = sprintf('%s_lfpPower.mat', fbase);
        ESIsave(fullfile(savefile, fname), 'lfpPower');
        clear lfpPower
        
        if allCfg.runErrorBars
            % Get event
            cfg = [];
            cfg.trials = trialsChosen;
            cfg.toilim = allCfg.timeSustained;
            lfpSel = ft_redefinetrial(cfg, lfpClean);
            
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
            
            fname = sprintf('%s_lfpPowerVar.mat', fbase);
            ESIsave(fullfile(savefile, fname), 'lfpPowerVar');
            clear TFRmultEvent
        end
    end

%% Connectivity
    function runConnectivity(lfpSel, muaSel)
        muaSel.label = strcat('muax-', muaSel.label);
        adata = ft_appenddata([], lfpSel, muaSel);
        
        %ppc auto
        cfg            = [];
        cfg.output     = 'powandcsd';
        cfg.method     = 'mtmfft';
        cfg.foi        = 10:2:120;
        cfg.taper      = 'dpss';
        cfg.tapsmofrq  = 7*ones(1,length(cfg.foi));
        cfg.keeptrials = 'yes';
        cfg.channel    = 'all';
        cfg.channelcmb = 'all';
        freq           = ft_freqanalysis(cfg, adata);
        
        cfg            = [];
        cfg.method     = 'coh';
        muaxCoherence  = ft_connectivityanalysis(cfg, freq);
        fname = sprintf('cond%02d_muaxCoherence.mat', cnd);
        ESIsave(fullfile(savefile, fname), 'muaxCoherence');
        clear muaxCoherence
    end

%% SFC
    function runSFC(spikeClean, lfpSel)
        cfg = [];
        cfg.trials = trialsClean(trialsChosen);
        cfg.toilim = allCfg.timeSustained;
        spikeSel = ft_spike_select(cfg, spikeClean);
        trl = unique([spikeSel.trial{:}]);
        spike.trialtime = repmat(allCfg.timeSustained, length(trl), 1);
        for ch = 1:nChan
            for tr = 1:length(trl)
                spikeSel.trial{ch}(spikeSel.trial{ch}==trl(tr)) = tr;
            end
        end
        spikeSel.label = strcat('label ',spikeSel.label);
        lfpSel.cfg.trl = repmat(allCfg.timeSustained, length(trl), 1);
        spikeSel.timestamp = spikeSel.time;
        
        % Then do it
        cfg = [];
        cfg.method = 'mtmconvol';
        cfg.foi    = 2:2:120;
        cfg.t_ftimwin = 7./cfg.foi; % watch out that the first frequency is not too long duration, so have to set to a minimum value
        twin = allCfg.timeSustained(2)-allCfg.timeSustained(1);
        cfg.t_ftimwin(cfg.t_ftimwin > twin) = twin;
        cfg.taper     = 'hanning';
        [sts] = ft_spiketriggeredspectrum(cfg, lfpSel, spikeSel); % for the event period
        
        % Get PPC
        cfg = [];
        cfg.method = 'ppc1'; % use ppc1
        for k = 1:nChan
            cfg.spikechannel = k;
            stSpec(k) = ft_spiketriggeredspectrum_stat(cfg, sts);
        end
        %         if ~allCfg.flagSecond
        fname = sprintf('%s_stSpec.mat', fbase);
        %         else
        %             fname = sprintf('%s_stSpec_Second.mat', fbase);
        %         end
        ESIsave(fullfile(savefile, fname), 'stSpec');
        clear stSpec
        
        if allCfg.runErrorBars
            for k = 1:nChan
                for rind = 1:sts.trial{k}(end)
                    cfg.spikechannel = k;
                    %         cfg.feedback = 'no';
                    cfg.trials = ~ismember([1:sts.trial{k}(end)],rind);
                    stSpecPerTrial(k, rind) = ft_spiketriggeredspectrum_stat(cfg, sts);
                end
            end
            fname = sprintf('cond%02d_stSpecPerTrial.mat', cnd);
            ESIsave(fullfile(savefile, fname), 'stSpecPerTrial');
            clear stSpecPerTrial
        end
    end

%% STA
    function runSTA(lfpSel, spikeSel)
        dataSel = ft_appendspike([], lfpSel, spikeSel);
        
        % Spike triggered average
        clear staPre
        cfg              = [];
        cfg.keeptrials   = 'no';
        cfg.latency      = allCfg.timeSustained;
        for k = 1:nChan
            cfg.spikechannel = spikeSel.label{k}; % first unit
            cfg.channel      = lfpSel.label(k); % first four chans
            cfg.timwin = [-0.2 0.2]; % take 400 ms
            cfg.feedback = 'no';
            sta(k) = ft_spiketriggeredaverage(cfg, dataSel);
        end
        fname = sprintf('cond%02d_sta.mat', cnd);
        ESIsave(fullfile(savefile, fname), 'sta');
        clear sta
    end
end
