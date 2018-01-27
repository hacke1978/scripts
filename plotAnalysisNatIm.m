clear all; close all; userpath('clear'); userpath('reset');
% Plot fieldtrip analysis for NatImFix paradigm
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Cem Uran
% cem.uran@esi-frankfurt.de
% 09/07/17
%%%%%%%%%%%%%%%%%%%%%%%%%%%
batchTag = '';  %analysis tag
% sesType = 'grating-ori'; % 'NatImFix' 'NatImSEQ' 'grating-ori' 'rfmapping-bar'
sesType = 'NatImFix';
dataDir = '/mnt/hpx/projects/MWNaturalPredict/Cem/Analysis';
saveDir = '/mnt/v7k/projects/MWNaturalPredict/Figures/';
% sessionList = {
% %     'hermes_20170613_fixation-naturalim_59',...
% %     'hermes_20170614_fixation-naturalim_61',...
% %     'hermes_20170616_fixation-naturalim_62',...
% %     'hermes_20170619_fixation-naturalim_63',...
% %     'hermes_20170620_fixation-naturalim_64_64b',...
% %     'hermes_20170623_fixation-naturalim_65',...
% %     'hermes_20170623_fixation-naturalim_66',...
% %     'hermes_20170626_fixation-naturalim_66',...
% %     'hermes_20170626_fixation-naturalim_67',...
% %     'hermes_20170627_fixation-naturalim_68_68b',...
% %     'hermes_20170628_fixation-naturalim_69',...
% %     'hermes_20170628_fixation-naturalim_70',...
% %     'hermes_20170629_fixation-naturalim_71',...
% %     'hermes_20170629_fixation-naturalim_72',...
% %     'hermes_20170630_fixation-naturalim_72',...
% %     'hermes_20170630_fixation-naturalim_73',...
% %     'hermes_20170703_fixation-naturalim_73',...
% %     'hermes_20171113_fixation-color-masksize_78_78b_78',...
% %     'hermes_20171115_fixation-color-masksize_79_79b',...
% %     'hermes_20171120_fixation-naturalim_80_80',...
% %     'hermes_20171120_fixation-naturalim_81',...
% %     'hermes_20171121_fixation-naturalim_82_82b',...
% %     'hermes_20171122_fixation-naturalim_83',...
% %     'hermes_20171123_fixation-naturalim_84_84b',...
% %     'hermes_20171124_fixation-naturalim_85_85',...
% %     'hermes_20171128_fixation-naturalim_88_88',...
% %     'hermes_20171128_fixation-naturalim_89_89_89_89',...
% %     'hermes_20171201_fixation-naturalim_90',...
% %     'hermes_20171204_fixation-naturalim_92',...
% %     'hermes_20171205_fixation-naturalim_93',...
% %     'hermes_20180110_fixation-naturalim_86',...
% %     'hermes_20180110_fixation-naturalim_96',...
%     };
% sessionList = {
%         'hermes_20171113_fixation-color-masksize_78_78b_78',...
%     'hermes_20171115_fixation-color-masksize_79_79b',...
% };
% sessionList = {
%     'hermes_20170425_fixation-grating-orientation-v2_1',...
%     'hermes_20171208_fixation-grating-orientation-v2_2'
% };
% sessionList = {
%     'hermes_20180110_fixation-naturalim_96',...
%     };

sessionList = {'hermes_20171127_fixation-naturalim_86'};
% sessionList = {'hermes_20171113_fixation-color-masksize_78_78b_78'};
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
        cfg{n}.do_trialPSTH = false;
        cfg{n}.do_lfpPower = false;
        cfg{n}.do_lfpPower2 = power2;
        cfg{n}.do_stSpec = true;
        cfg{n}.do_stSpecPerTrial = false;
        cfg{n}.withErrorBars = false;
        cfg{n}.onIm = false;
        cfg{n}.isOri = false;
        cfg{n}.isPred = false;
        cfg{n}.isCorr = false;
        cfg{n}.isGroup = false;
        cfg{n}.isOverlay = false;
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
