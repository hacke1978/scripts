function [trialEnd, trialStart, trialStimOn,  Cond, z] = findCompleteSEQNatIm(TDT)
if ~isempty(find(TDT.strobeCond == 255))
    TDT.tsDstr(TDT.strobeDstr == find(TDT.strobeCond == 255)) = [];
    TDT.strobeDstr(TDT.strobeDstr == find(TDT.strobeCond == 255)) = [];
end
TDT.tsCond(TDT.strobeCond==255)     = [];
TDT.strobeCond(TDT.strobeCond==255) = [];

TDT.tsEnd(TDT.strobeEnd==255)     = [];
TDT.strobeEnd(TDT.strobeEnd==255) = [];

TDT.tsEnd(TDT.strobeCond==0)      = [];
TDT.strobeEnd(TDT.strobeCond==0)  = [];
TDT.tsCond(TDT.strobeCond==0)     = [];
TDT.strobeCond(TDT.strobeCond==0) = [];

trialEnd     = TDT.tsEnd(TDT.strobeEnd~=1);
mistakes     = TDT.strobeEnd(TDT.strobeEnd~=1);

if length(TDT.tsCond)~=length(trialEnd)
    length(TDT.tsCond)
    length(trialEnd)
    error('start and end trial vectors have different sizes')
end

%
trialStart   = TDT.tsCond;
trialEnd     = trialEnd;
trialStimOn = repmat({nan(1, 2)}, length(trialEnd), 1);
stimOn = cellfun(@(x) TDT.tsDstr(find(x == TDT.strobeDstr)), ...
                 num2cell([1:length(TDT.strobeEnd)]+1), 'UniformOutput', false);
for ii=1:length(trialEnd)
    trialStimOn{ii}(:, 1:length(stimOn{ii})) = stimOn{ii};
end
trialStimOn = vertcat(trialStimOn{:});
Cond         = TDT.strobeCond;
z = mistakes;
