clear all; close all; userpath('clear'); userpath('reset');
% Run analysis for RFmapping paradigm
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Cem Uran
% cem.uran@esi-frankfurt.de
% 09/07/17
%%%%%%%%%%%%%%%%%%%%%%%%%%%
batchTag = '';
sesType = 'rfmapping-bar'; % 'NatImSEQ' 'grating-ori' 'rfmapping-bar'
monkeyName = 'Hermes'; % 'Hermes' 'Isis'
dirMonkey = fullfile('/mnt/hpx/projects/MWNaturalPredict', monkeyName, sesType);
saveMonkey = fullfile('/mnt/hpx/projects/MWNaturalPredict/Cem/Analysis', sesType, monkeyName);
% addpath('/mnt/hpx/slurm/uranc/fieldtrip/');

% sessionList = {...
% %   'ares014a04',...
%   'ares026a04',...
%   'ares030a03',...
% %     };
% sessionList = {
%     'hermes_20170808_rfmapping-bar_1'
% };
sessionList = {
    'hermes_20180110_rfmapping-bar_3',...
    'hermes_20180110_rfmapping-bar_4'};
saveMonkey = fullfile(saveMonkey, sessionList);
dirMonkey = fullfile(dirMonkey, sessionList);

%% Functions
cfg = [];
n = 0;
for ses = 1:length(sessionList)
    n = n+1;
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
    cfg{n}.print = true;
end

% Run analysis
tdt_RFanalysis_AP(cfg, 'local');