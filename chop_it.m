function [data] = chop_it(lfp, trialEnd, trialStart, Cond, taccept, preStim)

if length(preStim) == 1
    % Chop it - LFP
    data.fsample = lfp.data.fsample;
    data.cfg = lfp.data.cfg;
%     data.sampleinfo = lfp.data.sampleinfo;
    data.label = lfp.data.label;
    data.trialinfo = [[1:length(Cond)]' taccept Cond trialStart trialEnd]; % add in trialno
    time = lfp.data.time{1};
    for ii = 1:length(Cond)
        trialChosen = trialStart(ii) <= time & time <= trialEnd(ii);
        data.trial{ii} = lfp.data.trial{1}(:, trialChosen);
        data.time{ii} = time(trialChosen)-trialStart(ii)-preStim;
        data.sampleinfo{ii} = round(data.fsample*[trialStart(ii) trialEnd(ii) trialStart(ii)+preStim]);
    end
    data.sampleinfo = vertcat(data.sampleinfo{:});
else
    % Chop it - LFP
    data.fsample = lfp.data.fsample;
    data.cfg = lfp.data.cfg;
%     data.sampleinfo = lfp.data.sampleinfo;
    data.label = lfp.data.label;
    data.trialinfo = [[1:length(Cond)]' taccept Cond trialStart trialEnd]; % add in trialno
    time = lfp.data.time{1};
    for ii = 1:length(Cond)
        trialChosen = trialStart(ii) <= time & time <= trialEnd(ii);
        data.trial{ii} = lfp.data.trial{1}(:, trialChosen);
        data.time{ii} = time(trialChosen)-preStim(ii);
        data.sampleinfo{ii} = round(data.fsample*[trialStart(ii) trialEnd(ii) preStim(ii)]);
    end
    data.sampleinfo = vertcat(data.sampleinfo{:});
end