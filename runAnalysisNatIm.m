clear all; close all; userpath('clear'); userpath('reset');
% Run fieldtrip analysis for NatImFix paradigm
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Cem Uran
% cem.uran@esi-frankfurt.de
% 09/07/17
%%%%%%%%%%%%%%%%%%%%%%%%%%%
batchTag = '4sigma';  %analysis tag
sesType = 'NatImFix'; % 'NatImSEQ' 'grating-ori' 'rfmapping-bar'
monkeyName = 'Hermes'; % 'Hermes' 'Isis'
dirMonkey = fullfile('/mnt/hpx/projects/MWNaturalPredict', monkeyName, sesType);
saveMonkey = fullfile('/mnt/hpx/projects/MWNaturalPredict/Cem/Analysis', sesType, monkeyName);

% sessionList = dir(fullfile(dirMon, '*fixation-naturalim_*'));
% % sessionList = vertcat({sessionList.name});

sessionList = {...
    % 'ares022a03', ...
    % 'ares022a04', ...
    % 'ares023a02', ...
    % 'ares024a01', ...
    % 'ares025a03', ...
    % 'ares026a01', ...
    % 'ares026a02', ...
    % 'ares027a01', ...
    % 'ares027a02', ...
    % 'ares028a01', ...
    % 'ares029a01', ...
    % 'ares030a01', ...
    % 'ares030a02', ...
    % 'ares033a02', ...
    %     'ares033a03', ...
    %     'ares034a03', ...
    'hermes_20170808_rfmapping-bar_1'
    };

%% Functions
saveMonkey = fullfile(saveMonkey, sessionList);
dirMonkey = fullfile(dirMonkey, sessionList);

cfg = [];
n = 0;
for ses = 1:length(sessionList)
    n = n + 1;
    
    % get monkey dirs
    monkeyName = sessionList{ses}(1:4);
    if strcmp(monkeyName, 'herm'); monkeyName = 'Hermes'; end
    if strcmp(monkeyName, 'ares'); monkeyName = 'Ares'; end
    if strcmp(monkeyName, 'isis'); monkeyName = 'Isis'; end
    newSave = saveMonkey{ses};
    newDir = dirMonkey{ses};
    
    if ~exist(newSave, 'dir')
        mkdir(newSave); mkdir(fullfile(newSave, 'data'));
    end
    
    % set dirs
    cfg{n}.inputfile = newDir;
    cfg{n}.outputfile = newSave;
    cfg{n}.name = monkeyName;
    cfg{n}.type = sesType;
    cfg{n}.tag = batchTag;
    
    % set processing
    cfg{n}.filterLineNoise = ~strcmp(monkeyName, 'hermes');
    cfg{n}.getBadChannels = false;
    
    % set analysis options
    cfg{n}.runBaseline = false;
    cfg{n}.runErrorBars = false;
    cfg{n}.runPSTH = true;
    cfg{n}.runTimelockLFP = true;
    cfg{n}.runTimelockMUAX = true;
    cfg{n}.runConnectivity = false;
    cfg{n}.runTFR = false;
    cfg{n}.runSTA = false;
    cfg{n}.runSFC = false;
    
    % params
    cfg{n}.timeAll = [-0.5 3];
    cfg{n}.timeEvent = [0 3];
    cfg{n}.timeBaseline = [-0.5 0];
    cfg{n}.timeSustained = [0.25 0.75];
end

%% Run analysis
tdt_analysis_AP(cfg, 'local');
