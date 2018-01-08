
load(fullfile(filename, pathLFP.name), '-mat');
lfp(1).data = data; clear data;
load(fullfile(filename, pathMUAX.name), '-mat');
muax(1).data = data; clear data;
load(fullfile(filename, pathSpike.name), '-mat');
spike.timestamp = spike.time;
spikes(1).spike = spike;

extList = {
'hermes_20171211_fixation-grating-orientation-v2_3',...
'hermes_20171212_fixation-grating-orientation-v2_4',...
};

for ex=1:length(extList)
    lfp(ex+1) = load(fullfile(filename, pathLFP.name), '-mat');
    muax(ex+1) = load(fullfile(filename, pathMUAX.name), '-mat');
    spikes(ex+1) = load(fullfile(filename, pathSpike.name), '-mat');
    spikes(ex+1).spike.timestamp = spikes(ex+1).spike.time;
end

% % append two sesions
lfpAll = ft_appenddata([], lfp(1).data, lfp(2).data, lfp(3).data);
muaAll = ft_appenddata([], muax(1).data, muax(2).data, muax(3).data);
spikeAll = ft_appenddata([], spikes(1).spike, spikes(2).spike, spikes(3).spike);

spikeAll = spikes(1).spike;
for ii=2:3
spikeAll.trialinfo = [spikeAll.trialinfo; spikes(ii).spike.trialinfo];
spikeAll.sampleinfo = [spikeAll.sampleinfo; spikes(ii).spike.sampleinfo];
spikeAll.trialtime = [spikeAll.trialtime; spikes(ii).spike.trialtime];
spikeAll.cond = [spikeAll.cond; spikes(ii).spike.cond];
spikeAll.taccept = [spikeAll.taccept; spikes(ii).spike.taccept];
nTrials = max(spikeAll.trial{1});
for ch = 1:128
    spikeAll.time{ch} = [spikeAll.time{ch}; spikes(ii).spike.time{ch}];
    spikeAll.timestamp{ch} = [spikeAll.timestamp{ch}; spikes(ii).spike.timestamp{ch}];
    spikeAll.trial{ch} = [spikeAll.trial{ch}; spikes(ii).spike.trial{ch}+nTrials];
end
end
% clear mua lfp
% mua.data = muaAll; lfp.data = lfpAll;