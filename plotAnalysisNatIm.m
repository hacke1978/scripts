clear all; close all; userpath('clear'); userpath('reset');
% Plot fieldtrip analysis for NatImFix paradigm
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Cem Uran
% cem.uran@esi-frankfurt.de
% 09/07/17
%%%%%%%%%%%%%%%%%%%%%%%%%%%
batchTag = '';  %analysis tag
sesType = 'grating-ori'; % 'NatImSEQ' 'grating-ori' 'rfmapping-bar'

dataDir = '/mnt/hpx/projects/MWNaturalPredict/Cem/Analysis';
saveDir = '/mnt/v7k/projects/MWNaturalPredict/Figures/';
sessionList = {...
    %     'ares027a01', ...
    %     'ares027a02', ...
    %     'ares028a01', ...
    %     'ares029a01', ...
    %     'ares030a01', ...
    %     'ares030a02', ...
    %     'ares033a02',...
    %     'ares033a03',...
    %     'ares034a03',...
    % 'hermes_20171201_fixation-naturalim_90'
%     'ares034a03',...
%     'ares034a04',...
    'hermes_20171208_fixation-grating-orientation-v2_2'
    };
%% Functions

layType = 'stimuli';
powerLib = [false];
for power2 = powerLib
    cfg = [];
    n = 0;
    for ses = 1:length(sessionList)
        n = n+1;
        
        % get monkey dirs
        monkeyName = sessionList{ses}(1:4);
        if strcmp(monkeyName, 'herm'); monkeyName = 'Hermes'; end
        if strcmp(monkeyName, 'ares'); monkeyName = 'Ares'; end
        if strcmp(monkeyName, 'isis'); monkeyName = 'Isis'; end
        
        newDir = fullfile(dataDir, sesType, monkeyName, sessionList{ses});
        newSave = fullfile(saveDir, sesType, monkeyName, sessionList{ses});
        
        if ~exist(newSave, 'dir')
            mkdir(newSave); mkdir(fullfile(newSave, 'data'));
        end
        
        % set dirs
        cfg{n}.inputfile = newDir;
        cfg{n}.outputfile = newSave;
        cfg{n}.rfdir = fullfile(newDir);
        cfg{n}.name = monkeyName;
        cfg{n}.type = sesType;
        cfg{n}.tag = batchTag;
        
        % don't
        cfg{n}.do_muaxPPCSlepian = false;
        cfg{n}.do_muaxPPCHanning = false;
        cfg{n}.do_muaxCoherence = false;
        cfg{n}.do_sta = false;
        
        % do
        cfg{n}.do_timelockLFP = false;
        cfg{n}.do_timelockMUAX = false;
        cfg{n}.do_trialPSTH = true;
        cfg{n}.do_lfpPower = false;
        cfg{n}.do_lfpPower2 = power2;
        cfg{n}.do_stSpec = false;
        cfg{n}.do_stSpecPerTrial = false;
        cfg{n}.withErrorBars = false;
        cfg{n}.onIm = false;
        cfg{n}.isOri = false;
        cfg{n}.isPred = false;
        cfg{n}.isCorr = false;
        cfg{n}.isGroup = false;
        cfg{n}.window_rate = [0.05 0.15; 0.25 0.75];
        cfg{n}.window_freq = [40 70; 80 120];
        
        % set analysis options
        cfg{n}.corrLayout = 'channels';
        cfg{n}.layout = layType;
        cfg{n}.normalize = true;
        cfg{n}.print = true;
        cfg{n}.ext = 'png';
    end
    
    % Run analysis
    tdt_figures_AP(cfg, 'local');
end
