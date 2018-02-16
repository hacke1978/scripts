function [outfile] = runSFC_STA(allCfg, varargin)
if nargin > 1
    if strcmp(allCfg.type, 'grating-ori')
        [lfp, spike, taccept, caccept, Cond, CondSF] = deal(varargin{:});
    else
        [lfp, spike, taccept, caccept, Cond] = deal(varargin{:});
    end
    allCfg.spikeCh = find(caccept);
else
    if strcmp(allCfg.type, 'grating-ori')
        [lfp, ~, spike, taccept, caccept, Cond, CondSF] = load_preproc_data(allCfg);
    else
        [lfp, ~, spike, taccept, caccept, Cond] = load_preproc_data(allCfg);
    end
end
savefile = allCfg.outputfile;
outSFC = []; outSTA = []; outfile = [];
if allCfg.runErrorBars
    outSFCVar = [];
end
%% Should be the same at this point
% choose channels
% lfp.data.sampleinfo = [ lfp.data.sampleinfo lfp.data.sampleinfo(:, 1)];
cfg = [];
cfg.channel = find(caccept);
lfp = ft_selectdata(cfg, lfp.data);

cfg = [];
cfg.spikechannel = allCfg.spikeCh;
spike = ft_spike_select(cfg, spike);

% MOMENTARY FIX
lfp.hdr.Fs = lfp.fsample;
lfp.hdr.FirstTimeStamp = 0;
lfp.hdr.TimeStampPerSample = round(spike.hdr.Fs/lfp.fsample);
lfp.cfg.trl = lfp.sampleinfo;

spike.label = strcat(spike.label, '-spike');
data_all = ft_appendspike([], lfp, spike);

if strcmp(allCfg.type, 'grating-ori')
    condLib = unique(Cond);
    condLibSF = unique(CondSF);
    
    % then do stuff
    cfg              = [];
    cfg.keeptrials   = 'no';
    cfg.latency      = allCfg.timeSustained;
    for k = 1:length(allCfg.spikeCh)
        % SFC
        cfg = [];
        % cfg.method = 'mtmconvol';
        cfg.method = 'mtmfft';
        cfg.foi    = 2:2:120;
        cfg.t_ftimwin = 7./cfg.foi; % watch out that the first frequency is not too long duration, so have to set to a minimum value
        twin = allCfg.timeSustained(2)-allCfg.timeSustained(1);
        cfg.t_ftimwin(cfg.t_ftimwin>twin) = twin;
        cfg.taper     = 'hanning';
        [sts] = ft_spiketriggeredspectrum(cfg, data_all);
        
        for cnd = 1:length(condLib)
            for cnf = 1:length(condLibSF)
                
                %fname
                fbase = sprintf('cond%02d_%.2f_%s', condLib(cnd), condLibSF(cnf), allCfg.tag);
                
                % trials
                trialsChosen = Cond==condLib(cnd) & CondSF==condLibSF(cnf) & taccept==1;
                ch = allCfg.spikeCh(k);
                if allCfg.runSTA
                    cfg.trials = trialsChosen;
                    cfg.spikechannel = spike.label{k}; % first unit
                    cfg.channel      = lfp.label(k); % first four chans
                    cfg.timwin = [-0.2 0.2]; % take 400 ms
                    cfg.feedback = 'no';
                    stAv = ft_spiketriggeredaverage(cfg, data_all);
                    if ~allCfg.flagSecond
                        fname = sprintf('%sstAv_ch%02d.mat', fbase, ch);
                        outSTA = [outSTA; {fullfile(savefile, fname)}];
                        ESIsave(fullfile(savefile, fname), 'stAv');
                    else
                        stAv_Second = stAv;
                        fname = sprintf('%sstAv_Second_ch%02d.mat', fbase, ch);
                        outSTA = [outSTA; {fullfile(savefile, fname)}];
                        ESIsave(fullfile(savefile, fname), 'stAv_Second');
                    end
                    clear stAv
                end
                
                if allCfg.runSFC
                    cfg = [];
                    cfg.method = 'ppc1'; % use ppc1
                    cfg.trials = trialsChosen;
                    cfg.spikechannel = k;
                    stSpec = ft_spiketriggeredspectrum_stat(cfg, sts);
                    
                    % Save baseline
                    if ~allCfg.flagSecond
                        fname = sprintf('%sstSpec_ch%02d.mat', fbase, ch);
                        outSFC = [outSFC; {fullfile(savefile, fname)}];
                        ESIsave(fullfile(savefile, fname), 'stSpec');
                        clear stSpec
                    else
                        stSpec_Second = stSpec;
                        fname = sprintf('%sstSpec_Second_ch%02d.mat', fbase, ch);
                        outSFC = [outSFC; {fullfile(savefile, fname)}];
                        ESIsave(fullfile(savefile, fname), 'stSpec_Second');
                        clear stSpec_Second
                    end
                    if allCfg.runErrorBars
                        cfg = [];
                        cfg.method = 'ppc1'; % use ppc1
                        cfg.spikechannel = spike.label{k}; % first unit
                        trAll = find(trialsChosen);
                        for rind = 1:length(trAll)
                            tind = trAll(rind);
                            cfg.trials = trAll(~ismember([trAll], tind));
                            stSpecPerTrial(rind) = ft_spiketriggeredspectrum_stat(cfg, sts);
                        end
                        if ~allCfg.flagSecond
                            fname = sprintf('%sstSpecPerTrial_ch%02d.mat', fbase, ch);
                            outSFCVar = [outSFCVar; {fullfile(savefile, fname)}];
                            ESIsave(fullfile(savefile, fname), 'stSpecPerTrial');
                            clear stSpecPerTrial
                        else
                            stSpecPerTrial_Second = stSpecPerTrial;
                            fname = sprintf('%sstSpecPerTrial_Second_ch%02d.mat', fbase, ch);
                            outSFCVar = [outSFCVar; {fullfile(savefile, fname)}];
                            ESIsave(fullfile(savefile, fname), 'stSpecPerTrial_Second');
                            clear stSpecPerTrial_Second
                        end
                    end
                end
            end
        end
    end
else
    condLib = unique(Cond);
    % then do stuff
    cfg              = [];
    cfg.keeptrials   = 'no';
    cfg.latency      = allCfg.timeSustained;
    for k = 1:length(allCfg.spikeCh)
        for cnd = 1:length(condLib)
            %fname
            fbase = sprintf('cond%02d_%s', condLib(cnd), allCfg.tag);
            % trials
            trialsChosen = (Cond==condLib(cnd)) & taccept==1;
            ch = allCfg.spikeCh(k);
            if allCfg.runSTA
                cfg.trials = trialsChosen;
                cfg.spikechannel = spike.label{k}; % first unit
                cfg.channel      = lfp.label; % first four chans
                cfg.timwin = [-0.2 0.2]; % take 400 ms
                cfg.feedback = 'no';
                stAv = ft_spiketriggeredaverage(cfg, data_all);
                if ~allCfg.flagSecond
                    fname = sprintf('%sstAv_ch%02d.mat', fbase, ch);
                    outSTA = [outSTA; {fullfile(savefile, fname)}];
                    ESIsave(fullfile(savefile, fname), 'stAv');
                else
                    stAv_Second = stAv;
                    fname = sprintf('%sstAv_Second_ch%02d.mat', fbase, ch);
                    outSTA = [outSTA; {fullfile(savefile, fname)}];
                    ESIsave(fullfile(savefile, fname), 'stAv_Second');
                end
                clear sta
            end
            if allCfg.runSFC
                % SFC
                cfg = [];
                % cfg.method = 'mtmconvol';
                cfg.spikechannel = spike.label{k}; % first unit
                cfg.method = 'mtmfft';
                cfg.foi    = 2:2:120;
                cfg.t_ftimwin = 7./cfg.foi; % watch out that the first frequency is not too long duration, so have to set to a minimum value
                twin = allCfg.timeSustained(2)-allCfg.timeSustained(1);
                cfg.t_ftimwin(cfg.t_ftimwin>twin) = twin;
                cfg.taper     = 'hanning';
                [sts] = ft_spiketriggeredspectrum(cfg, data_all);
                
                cfg = [];
                cfg.method = 'ppc1'; % use ppc1
                cfg.trials = trialsChosen;
                cfg.spikechannel = spike.label{k}; % first unit
                stSpec = ft_spiketriggeredspectrum_stat(cfg, sts);
                
                % Save baseline
                if ~allCfg.flagSecond
                    fname = sprintf('%sstSpec_ch%02d.mat', fbase, ch);
                    outSFC = [outSFC; {fullfile(savefile, fname)}];
                    ESIsave(fullfile(savefile, fname), 'stSpec');
                    clear stSpec
                else
                    stSpec_Second = stSpec;
                    fname = sprintf('%sstSpec_Second_ch%02d.mat', fbase, ch);
                    outSFC = [outSFC; {fullfile(savefile, fname)}];
                    ESIsave(fullfile(savefile, fname), 'stSpec_Second');
                    clear stSpec_Second
                end
                
                if allCfg.runErrorBars
                    cfg = [];
                    cfg.method = 'ppc1'; % use ppc1
                    cfg.spikechannel = spike.label{k}; % first unit
                    trAll = find(trialsChosen);
                    for rind = 1:length(trAll)
                        tind = trAll(rind);
                        cfg.trials = trAll(~ismember([trAll], tind));
                        stSpecPerTrial(rind) = ft_spiketriggeredspectrum_stat(cfg, sts);
                    end
                    if ~allCfg.flagSecond
                        fname = sprintf('%sstSpecPerTrial_ch%02d.mat', fbase, ch);
                        outSFCVar = [outSFCVar; {fullfile(savefile, fname)}];
                        ESIsave(fullfile(savefile, fname), 'stSpecPerTrial');
                        clear stSpecPerTrial
                    else
                        stSpecPerTrial_Second = stSpecPerTrial;
                        fname = sprintf('%sstSpecPerTrial_Second_ch%02d.mat', fbase, ch);
                        outSFCVar = [outSFCVar; {fullfile(savefile, fname)}];
                        ESIsave(fullfile(savefile, fname), 'stSpecPerTrial_Second');
                        clear stSpecPerTrial_Second
                    end
                end
            end
        end
    end
end

%% get the filenames
outfile.staFiles = outSTA;
outfile.sfcFiles = outSFC;
if allCfg.runErrorBars
    outfile.sfcFilesVar = outSFCVar;
end