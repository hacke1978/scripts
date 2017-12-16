function [varargout] = getSnipsForSession(sesName, sesType)
addpath('/mnt/hpx/slurm/uranc/TDT2MAT')
addpath('/mnt/hpx/slurm/uranc/TDT2MAT/TDT_MichaelStephan_Add2Path/')

% Gets the trial info to make snippets
global outputDirSes
TDTFilename = sprintf('%s/Experiments_%s', outputDirSes, sesName);
% MATFilename = sprintf('/mnt/hpx/projects/MWNaturalPredict/Ares/%s/%s', sesName, sesName);

disp('load TDT')
DataTDT    = ft_TDT(TDTFilename);
nrChannels = 32;

TDT.starttime = DataTDT.starttime;
for iStrobes = 1: length(DataTDT.strobeStores(:))
    if strcmp(DataTDT.strobeStores(1,iStrobes).StoreName,'Cond')
        TDT.strobeCond  = DataTDT.strobeStores(1,iStrobes).strobes;
        TDT.tsCond      = DataTDT.strobeStores(1,iStrobes).timeStamps - DataTDT.starttime;
    elseif strcmp(DataTDT.strobeStores(1,iStrobes).StoreName,'Trgt')
        TDT.strobeTrgt  = DataTDT.strobeStores(1,iStrobes).strobes;
        TDT.tsTrgt      = DataTDT.strobeStores(1,iStrobes).timeStamps - DataTDT.starttime;
    elseif strcmp(DataTDT.strobeStores(1,iStrobes).StoreName,'Tend')
        TDT.strobeEnd   = DataTDT.strobeStores(1,iStrobes).strobes;
        TDT.tsEnd       = DataTDT.strobeStores(1,iStrobes).timeStamps - DataTDT.starttime;
        %elseif strcmp(DataTDT.strobeStores(1,iStrobes).StoreName,'Tnum')
        %TDT.strobeTrlNo = DataTDT.strobeStores(1,iStrobes).strobes;
    elseif strcmp(DataTDT.strobeStores(1,iStrobes).StoreName,'Dstr')
        TDT.strobeDstr  = DataTDT.strobeStores(1,iStrobes).strobes;
        TDT.tsDstr      = DataTDT.strobeStores(1,iStrobes).timeStamps - DataTDT.starttime;
    elseif strcmp(DataTDT.strobeStores(1,iStrobes).StoreName,'Chng')
        TDT.strobeChng  = DataTDT.strobeStores(1,iStrobes).strobes;
        TDT.tsChng      = DataTDT.strobeStores(1,iStrobes).timeStamps - DataTDT.starttime;
    end
end

% Get trial times
if strcmp(sesType, 'NatImFix')
    [trialEnd, trialStart, trialStimOn, Cond, taccept] = findCompleteFlashVinck(TDT);
    varargout{1} = trialEnd;
    varargout{2} = trialStart;
    varargout{3} = trialStimOn;
    varargout{4} = Cond;
    varargout{5} = taccept;
elseif strcmp(sesType, 'NatImSEQ')
    if strcmp(sesName, 'ares015a02')
        TDT.strobeCond = [255; TDT.strobeCond];
        TDT.tsCond = [0; TDT.tsCond];
    end
    [trialEnd, trialStart, trialStimOn, Cond, taccept] = findCompleteSEQNatIm(TDT);
    varargout{1} = trialEnd;
    varargout{2} = trialStart;
    varargout{3} = trialStimOn;
    varargout{4} = Cond;
    varargout{5} = taccept;
elseif strcmp(sesType, 'rfmapping-bar')
    [trialEnd, trialStart, trialStimOn, Cond, taccept] = findCompleteFlashVinck(TDT);
    varargout{1} = trialEnd;
    varargout{2} = trialStart;
    varargout{3} = trialStimOn-0.5;
    varargout{4} = Cond;
    varargout{5} = taccept;
else
end
