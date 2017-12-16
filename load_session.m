function [allFiles] = load_session(allCfg, filename)

% look for files
path_lfpPower_Second = dir(fullfile(filename, allCfg.tag, '*lfpPower_Second.mat'));
path_lfpPower = dir(fullfile(filename, allCfg.tag, '*_lfpPower.mat'));
path_lfpPower2 = dir(fullfile(filename, allCfg.tag, '*_lfpPower2.mat'));
path_lfpPowerVar = dir(fullfile(filename, allCfg.tag, '*_lfpPowerVar.mat'));
path_muaxPPCSlepian = dir(fullfile(filename, allCfg.tag, '*_muaxPPCSlepian.mat'));
path_muaxPPCHanning = dir(fullfile(filename, allCfg.tag, '*_muaxPPCHanning.mat'));
path_muaxCoherence = dir(fullfile(filename, allCfg.tag, '*_muaxCoherence.mat'));
path_sta = dir(fullfile(filename, allCfg.tag, '*_sta.mat'));
path_stSpec = dir(fullfile(filename, allCfg.tag, '*_stSpec.mat'));
path_stSpecPerTrial = dir(fullfile(filename, allCfg.tag, '*_stSpecPerTrial.mat'));
path_timelockLFP = dir(fullfile(filename, allCfg.tag, '*_timelockLFP.mat'));
path_timelockMUAX = dir(fullfile(filename, allCfg.tag, '*_timelockMUAX.mat'));
path_trialPSTH = dir(fullfile(filename, allCfg.tag, '*_trialPSTH.mat'));
path_trialPSTHVar = dir(fullfile(filename, allCfg.tag, '*_trialPSTHVar.mat'));
path_powerBaseline = dir(fullfile(filename, allCfg.tag, '*_powerBaseline.mat'));
path_stSpecBaseline = dir(fullfile(filename, allCfg.tag, '*_stSpecBaseline*'));

% also the baseline
allCfg.path_powerBaseline = ~isempty(path_powerBaseline);
allCfg.path_stSpecBaseline = ~isempty(path_stSpecBaseline);

% prepare the config
allCfg.path_lfpPower_Second = ~isempty(path_lfpPower_Second) && allCfg.do_lfpPower_Second;
allCfg.path_lfpPower = ~isempty(path_lfpPower) && allCfg.do_lfpPower;
allCfg.path_lfpPower2 = ~isempty(path_lfpPower2) && allCfg.do_lfpPower2;
allCfg.path_lfpPowerVar = ~isempty(path_lfpPowerVar) && allCfg.isShaded;
allCfg.path_muaxPPCSlepian = ~isempty(path_muaxPPCSlepian) && allCfg.do_muaxPPCSlepian;
allCfg.path_muaxPPCHanning = ~isempty(path_muaxPPCHanning) && allCfg.do_muaxPPCHanning;
allCfg.path_muaxCoherence = ~isempty(path_muaxCoherence) && allCfg.do_muaxCoherence;
allCfg.path_timelockLFP = ~isempty(path_timelockLFP) && allCfg.do_timelockLFP;
allCfg.path_timelockMUAX = ~isempty(path_timelockMUAX) && allCfg.do_timelockMUAX;
allCfg.path_trialPSTH = ~isempty(path_trialPSTH) && allCfg.do_trialPSTH;
allCfg.path_trialPSTHVar = ~isempty(path_trialPSTHVar) && allCfg.isShaded;
allCfg.path_sta = ~isempty(path_sta) && allCfg.do_sta;
allCfg.path_stSpec = ~isempty(path_stSpec) && allCfg.do_stSpec;
allCfg.path_stSpecPerTrial = ~isempty(path_stSpecPerTrial) && allCfg.do_stSpecPerTrial;

% check the length
fields = fieldnames(allCfg);
pind = find(strncmp(fields, 'path_', 4));
allLength = []; ind2del = [];
for ii = 1:length(pind)
    if allCfg.(fields{pind(ii)})
        fname = eval('fields{pind(ii)}');
        if strcmp(fname(end-3:end), 'line')
            assert(length(eval(fname)) == 1, 'Too many baseline files!')
        else
            allLength = [allLength  length(eval(fname))];
        end
    else
        ind2del = [ind2del ii];
    end
end
pind(ind2del) = [];
assert(length(unique(allLength)) == 1, 'Some files are missing!')
nFiles = length(pind);
nCond = unique(allLength);
% load all files
allFiles = [];
for nf = 1:nFiles
    for cnd = 1:nCond
        vname = (fields{pind(nf)});
        thisVar = eval(vname);
        if strcmp(vname(end-3:end), 'line')
            ESIload(fullfile(filename, thisVar.name));
            allFiles.(vname(6:end)) = eval(vname(6:end));
            break
        else
            ESIload(fullfile(filename, thisVar(cnd).name));
            allFiles(cnd).(vname(6:end)) = eval(vname(6:end));
        end
    end
end