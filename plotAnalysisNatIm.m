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
% sesType = 'NatImSEQ';
dataDir = '/mnt/hpx/projects/MWNaturalPredict/Cem/Analysis';
saveDir = '/mnt/v7k/projects/MWNaturalPredict/Figures/';
% sessionList = {...
%     'hermes_20170531_fixation-naturalim_49',...
%     'hermes_20170607_fixation-naturalim_57',...
%     'hermes_20170608_fixation-naturalim_58',...
%     'hermes_20170608_fixation-naturalim_58b',...
%     'hermes_20170612_fixation-naturalim_57',...
%     'hermes_20170613_fixation-naturalim_59',...
%     'hermes_20170614_fixation-naturalim_61',...
%     'hermes_20170616_fixation-naturalim_62',...
%     'hermes_20170619_fixation-naturalim_63',...
%     'hermes_20170620_fixation-naturalim_64',...
%     'hermes_20170620_fixation-naturalim_64b',...
%     'hermes_20170623_fixation-naturalim_65',...
%     'hermes_20170623_fixation-naturalim_66',...
%     'hermes_20170626_fixation-naturalim_66',...
%     'hermes_20170626_fixation-naturalim_67',...
%     'hermes_20170627_fixation-naturalim_68',...
%     'hermes_20170627_fixation-naturalim_68b',...
%     'hermes_20170628_fixation-naturalim_69',...
%     'hermes_20170628_fixation-naturalim_70',...
%     'hermes_20170629_fixation-naturalim_71',...
%     'hermes_20170629_fixation-naturalim_72',...
%     'hermes_20170630_fixation-naturalim_72',...
%     'hermes_20170630_fixation-naturalim_73',...
%     'hermes_20170703_fixation-naturalim_73',...
%     'hermes_20171120_fixation-naturalim_80',...
%     'hermes_20171120_fixation-naturalim_81',...
%     'hermes_20171121_fixation-naturalim_82',...
%     'hermes_20171121_fixation-naturalim_82b',...
%     'hermes_20171122_fixation-naturalim_83',...
%     'hermes_20171123_fixation-naturalim_84',...
%     'hermes_20171123_fixation-naturalim_84b',...
%     'hermes_20171124_fixation-naturalim_85',...
%     'hermes_20171127_fixation-naturalim_80',...
%     'hermes_20171128_fixation-naturalim_88',...
%     'hermes_20171128_fixation-naturalim_89',...
%     'hermes_20171129_fixation-naturalim_89',...
%     'hermes_20171130_fixation-naturalim_85',...
%     'hermes_20171130_fixation-naturalim_89',...
%     'hermes_20171201_fixation-naturalim_90',...
%     'hermes_20171204_fixation-naturalim_92',...
%     'hermes_20171204_fixation-naturalim_93',...
%     'hermes_20171205_fixation-naturalim_93',...
%     'hermes_20171206_fixation-naturalim_88',...
%     'hermes_20171206_fixation-naturalim_89',...
%     'hermes_20171207_fixation-naturalim_96',...
%     'hermes_20180110_fixation-naturalim_86',...
%     'hermes_20180110_fixation-naturalim_96',...
%     'hermes_20180118_fixation-naturalim_86B',...
%     'hermes_20180118_fixation-naturalim_86G',...
%     'hermes_20180118_fixation-naturalim_86R',...
%     'hermes_20180122_fixation-naturalim_86B',...
%     'hermes_20180122_fixation-naturalim_86G',...
%     'hermes_20180122_fixation-naturalim_86k',...
%     'hermes_20180122_fixation-naturalim_86R',...
%     'hermes_20180124_fixation-naturalim_96b',...
%     'hermes_20180123_fixation-naturalim_96',...
%     'hermes_20180124_fixation-naturalim_96',...
%     'hermes_20180125_fixation-naturalim_97',...
%     'hermes_20180125_fixation-naturalim_97b',...
%     'hermes_20180125_fixation-naturalim_98',...
%     'hermes_20180126_fixation-naturalim_96',...
%     'hermes_20180129_fixation-naturalim_101',...
%     'hermes_20180129_fixation-naturalim_102',...
%     'hermes_20180130_fixation-naturalim_102',...
%     'hermes_20180131_fixation-naturalim_102',...
%     'hermes_20180131_fixation-naturalim_103B',...
%     'hermes_20180131_fixation-naturalim_103Y',...
%         'ares025a03',...
%     'ares026a01',...
%     'ares026a02',...
%     'ares027a01',...
%     'ares027a02',...
%     'ares028a01',...
%     'ares029a01',...
%     'ares030a01',...
%     'ares030a02',...
%     'ares034a03',...
%     'ares034a04',...
%     'ares037a01',...
%     'ares038a01',...
%     'ares038a02',...
%     'ares039a01',...
%     'ares039a02',...
%     'ares040a01',...
%     'ares040a02',...
%     'ares040a03',...
%     'ares042a01',...
%     'ares042a02',...  
%     };
sessionList = {
%     'hermes_20170628_fixation-naturalim_69',...
    'hermes_20170607_fixation-naturalim_57',...
%     'hermes_20171113_fixation-color-masksize_78_78b_78',...
%     'hermes_20171115_fixation-color-masksize_79_79b',...
};
% sessionList = {
%     'hermes_20170425_fixation-grating-orientation-v2_1',...
% %     'hermes_20171208_fixation-grating-orientation-v2_2'
% };
% sessionList = {
% 'hermes_20180201_fixation-naturalim_103E',...
% 'hermes_20180201_fixation-naturalim_103K',...
% 'hermes_20180202_fixation-naturalim_103E',...
% 'hermes_20180202_fixation-naturalim_103G',...
% 'hermes_20180202_fixation-naturalim_103R',...
% };
% sessionList = {
%     'hermes_20170712_fixation-naturalim_60',...
%     'hermes_20170724_fixation-naturalim_70',...
%     'hermes_20170725_fixation-naturalim_71',...
%     'hermes_20180205_fixation-naturalim_102',...
%     'hermes_20180207_fixation-naturalim_103E',...
%     'hermes_20180207_fixation-naturalim_103G3',...
%     'hermes_20180207_fixation-naturalim_103K',...
%     'hermes_20180207_fixation-naturalim_103W'
%     };
% sessionList = {
%         'hermes_20180118_fixation-naturalim_86B',...
%     'hermes_20180118_fixation-naturalim_86G',...
%     'hermes_20180118_fixation-naturalim_86R',...
%     'hermes_20180122_fixation-naturalim_86B',...
%     'hermes_20180122_fixation-naturalim_86G',...
%     'hermes_20180122_fixation-naturalim_86k',...
%     'hermes_20180122_fixation-naturalim_86R',...
%     'hermes_20170628_fixation-naturalim_70',...
%     'hermes_20170630_fixation-naturalim_73',...
%     'hermes_20170703_fixation-naturalim_73',...
%     'hermes_20171127_fixation-naturalim_86',...
%     'hermes_20171127_fixation-naturalim_87',...
% };

% sessionList = {
%     'hermes_20170607_fixation-naturalim_57',...
%     'ares037a01',...
%     'ares038a01',...
%     'ares038a02',...
%     'ares033a02',...
%     'ares033a03',...
%     'ares039a01',...
%     'ares039a02',...
%     };
%
% sessionList = {
% %     'hermes_20171113_fixation-color-masksize_78',...
%     'hermes_20171114_fixation-color-masksize_78b',...
% %     'hermes_20171116_fixation-color-masksize_78',...
%     };
% sessionList = {
%     'hermes_20180129_fixation-naturalim_101',...
%     'hermes_20180110_fixation-naturalim_86',...
%     'hermes_20180129_fixation-naturalim_102',...
%     'hermes_20180126_fixation-naturalim_96'
% };
% sessionList = {'hermes_20180125_fixation-naturalim_97_97b',...
%                'hermes_20180110_fixation-naturalim_96'};
% sessionList = {'hermes_20171127_fixation-naturalim_86'};
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
        cfg{n}.do_sta = false;
        
        % do
        cfg{n}.gammaFitRange = [15:140];% 'all';
        cfg{n}.gammaFitType = 'exp'; % 'exp' 'poly' 'linear'
        %  cfg{n}.gammaPeak = 'all';
        %  cfg{n}.gammaPeak = [20:120];
        %  cfg{n}.gammaPeak = [35:100]; % example 35:80
        cfg{n}.do_muaxCoherence = false;
        cfg{n}.do_timelockLFP = false;
        cfg{n}.do_timelockMUAX = false;
        cfg{n}.do_trialPSTH = false;
        cfg{n}.do_lfpPower = true;
        cfg{n}.do_lfpPower2 = power2;
        cfg{n}.do_stSpec = false;
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
%     try
    tdt_figures_AP(cfg, 'local');
%     catch
%     end
end
