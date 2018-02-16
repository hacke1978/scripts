function analyze_session(allCfg)
addpath('/mnt/hpx/slurm/uranc/fieldtrip/');
% addpath('/mnt/v7k/projects/MWNaturalPredict/fieldtrip');
ft_defaults

if ~isfield(allCfg, 'filterLineNoise');    allCfg.filterLineNoise = false; end
if ~isfield(allCfg, 'getBadChannels');     allCfg.getBadChannels = false;  end
if ~isfield(allCfg, 'runPSTH');            allCfg.timelockLfp = false;     end
if ~isfield(allCfg, 'runTimelockLFP');     allCfg.timelockLfp = false;     end
if ~isfield(allCfg, 'runTimelockMUAX');    allCfg.timelockMuax = false;    end
if ~isfield(allCfg, 'runTFR');             allCfg.tfr = false;             end
if ~isfield(allCfg, 'save');               allCfg.save = false;            end


% check filename
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

% Load data - LFP, MUA, Spike
if strcmp(allCfg.type, 'grating-ori')
    [lfp, muax, spike, taccept, caccept, Cond, CondSF] = load_preproc_data(allCfg);
else
    [lfp, muax, spike, taccept, caccept, Cond] = load_preproc_data(allCfg);
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
CondAll = Cond;
Cond = Cond(taccept);
trialsClean = find(taccept);
condLib = unique(Cond);

if strcmp(allCfg.type, 'grating-ori')
    CondSFAll = CondSF;
    CondSF = CondSF(taccept);
    condLibSF = unique(CondSF);
    for sf=1:length(condLibSF)
        for cnd=1:length(condLib)
            cnd
            trialsChosen = find(Cond == condLib(cnd) & CondSF == condLibSF(sf));
            analyze_condition(allCfg, lfpClean, muaClean, spikeClean, trialsClean, trialsChosen, condLib(cnd), condLibSF(sf));
        end
    end
elseif strcmp(allCfg.type, 'NatImSEQ')
    for cnd=1:length(condLib)
        
        % choose trials
        trialsChosen = find(Cond==condLib(cnd));
        
        % first stim
        allCfg.flagSecond = false;
        allCfg.timeSustained = allCfg.timeSustained_First;
        analyze_condition(allCfg, lfpClean, muaClean, spikeClean, trialsClean, trialsChosen, condLib(cnd));
        
        % second stim
        allCfg.flagSecond = true;
        allCfg.timeSustained = allCfg.timeSustained_Second;
        analyze_condition(allCfg, lfpClean, muaClean, spikeClean, trialsClean, trialsChosen, condLib(cnd));
        
    end
else
    for cnd=1:length(condLib)
        trialsChosen = find(Cond==condLib(cnd))';
        length(trialsChosen)
        %             tok = strsplit(filename, '/');
        %             if strcmp(tok{end}, 'hermes_20170428_fixation-naturalim_15')
        %                 trialsChosen = trialsChosen(1:18);
        %             end
        analyze_condition(allCfg, lfpClean, muaClean, spikeClean, trialsClean, trialsChosen, condLib(cnd))
    end
end

% Save Baseline TFR
if (allCfg.runTFR && allCfg.runBaseline); runTFR_Baseline(allCfg, lfpClean); end

% send the rest to slurm switch calcLocation
if (allCfg.runSTA || allCfg.runSFC)
    % get the channels and organize
    cSel = find(caccept);
    calcLocation = 'slurm';
    spikeCfg = {};
    for ii=1:length(cSel)
        spikeCfg{ii} = allCfg;
        spikeCfg{ii}.spikeCh = cSel(ii);
    end
    
    switch calcLocation
        case 'slurm'
            license('inuse')
            out = slurmfun(@runSFC_STA, spikeCfg, 'stopOnError', false, 'partition', '16GB', 'useUserPath', true);
%             out = slurmfun(@runSFC_STA, spikeCfg, 'stopOnError', false, 'partition', '24GBL', 'useUserPath', true);
%             out = slurmfun(@analyze_session, allCfg, 'stopOnError', false, 'partition', '48GBL', 'useUserPath', true);
        case 'local'
            if strcmp(allCfg.type, 'grating-ori')
                out = cellfun(@runSFC_STA, {allCfg}, {lfp}, {spike}, {taccept}, {caccept}, {CondAll}, {CondSFAll});
            else
                out = cellfun(@runSFC_STA, {allCfg}, {lfp}, {spike}, {taccept}, {caccept}, {CondAll});
            end
    end
    merge_analysis(out, true)
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Subfunctions %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% TFR Baseline
function runTFR_Baseline(allCfg, lfpClean)
savename = allCfg.outputfile;
cfg = [];
cfg.trials = 'all';
cfg.toilim = allCfg.timeBaseline;
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
return
