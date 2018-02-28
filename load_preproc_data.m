function [varargout]  =  load_preproc_data(allCfg)

filename = allCfg.inputfile;
if ~iscell(filename); filename = {filename}; end;
% % merge sessions if any
% if iscell(allCfg.outputfile)
%     if strcmp(allCfg.name, 'Hermes')
%         tok = cellfun(@(x) strsplit(x, '_'), allCfg.outputfile, 'UniformOutput', false);
%         tok = vertcat(tok{:});
%         mergeName  = allCfg.outputfile{1};
%         for ii=2:length(allCfg.outputfile)
%             mergeName = [mergeName , sprintf('_%s', tok{ii, end})];
%         end
%     else
%         tok = cellfun(@(x) strsplit(x, '/'), allCfg.outputfile, 'UniformOutput', false);
%         tok = vertcat(tok{:});
%         mergeName  = allCfg.outputfile{1};
%         for ii=2:length(allCfg.outputfile)
%             mergeName = [mergeName , sprintf('_%s', tok{ii, end}(end-2:end))];
%         end
%     end
%     allCfg.outputfile = mergeName;
% end

% saveDir
if ~exist(allCfg.outputfile, 'dir')
    mkdir(allCfg.outputfile);
end

% Load the data - handles session merging
for ii=1:length(filename);
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
        nTrials = max(cellfun(@max, spike.trial));
        for ch = 1:length(spike.label)
            spike.time{ch} = [spike.time{ch}; spikeAdd.spike.time{ch}];
            spike.timestamp{ch} = [spike.timestamp{ch}; spikeAdd.spike.timestamp{ch}];
            spike.trial{ch} = [spike.trial{ch}; spikeAdd.spike.trial{ch}+nTrials];
        end
    end
end

if length(filename)>1
    
    % save lfp
    data = lfp.data;
    save([allCfg.outputfile '_chopped.lfp'],'data', '-v7.3');
    clear data;
    
    % save muax
    data = muax.data;
    save([allCfg.outputfile '_chopped.muax'],'data', '-v7.3');
    clear data;
    
    % save spike
    save([allCfg.outputfile '_chopped.spike'],'spike', '-v7.3');
end

% get conditions / channels
if strcmp(allCfg.name, 'Hermes')
    if strcmp(allCfg.type, 'grating-ori')
        %         taccept = (muax.data.trialinfo(:, 4)==1);
        taccept = (muax.data.trialinfo(:, 2)==1);
        CondSF = muax.data.trialinfo(:, 4);
    elseif strcmp(allCfg.type, 'rfmapping-bar')
        taccept = (true([size(muax.data.trialinfo, 1) 1]));
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

% varargout
varargout = [];
varargout{1} = lfp;
varargout{2} = muax;
varargout{3} = spike;
varargout{4} = taccept;
varargout{5} = caccept;
varargout{6} = Cond;
if strcmp(allCfg.type, 'grating-ori'); varargout{7} = CondSF; end
