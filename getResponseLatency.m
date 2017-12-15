function [gMap, gLatency] = getResponseLatency(psth, edges)
orient = 0:22.5:360-22.5;
respLatency = 0.01:0.001:0.1;
stimDuration = 3;
thisPSTH = vertcat(psth{:});
thisPSTH = [thisPSTH(9:end, :); thisPSTH(1:8, :)];
allMap = []; allMax = [];
for rind =1:length(respLatency)
    rDel = respLatency(rind);
    tsel = rDel < edges & edges < rDel+stimDuration;
    delPSTH = thisPSTH(:, tsel);
    allMap{rind} = back_project(delPSTH, orient);
    allMax = [allMax max(allMap{rind}(:))];
end
[~, ind] = max(allMax);
gMap = allMap{ind}; gLatency = respLatency(ind);