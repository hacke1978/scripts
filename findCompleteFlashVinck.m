function [trialEnd, trialStart, trialStimOn, Cond, z] = findCompleteFlashVinck(TDT)
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
%% not sure if this part is necessary
if length(TDT.tsCond)~=length(trialEnd)
    length(TDT.tsCond)
    length(trialEnd)
    error('start and end trial vectors have different sizes')
end

trialStart   = TDT.tsCond;
trialStimOn = trialStart + 0.5;
trialEnd     = trialEnd;
Cond         = TDT.strobeCond;
z = mistakes;