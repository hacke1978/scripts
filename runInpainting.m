% Run inpainting for
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Cem Uran
% cem.uran@esi-frankfurt.de
% 09/07/17
%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all; userpath('clear'); userpath('reset');
addpath('/mnt/hpx/opt/ESIsoftware/matlab/')
stimSetAll = '/mnt/v7k/home/uranc/workspace/VinckLab/Dataset/Processed/_6flickr/stimsetChosen';
stimDir = '/mnt/v7k/home/uranc/workspace/VinckLab/Dataset/Processed/_6flickr/stimset01';
setDir = '/mnt/v7k/projects/MWNaturalPredict/Dataset/YahooFlickr/All';
load(fullfile(stimSetAll, 'imSetHighSinger.mat'));
load(fullfile(stimSetAll, 'imSetLowSinger.mat'));
%%
name = 'ares';
saveDir = '/mnt/hpx/projects/MWNaturalPredict/Cem/Analysis/Inpainting/';
if strcmp(name, 'Hermes')
    load(fullfile('/mnt/hpx/projects/MWNaturalPredict/Cem/Analysis/RFmaps/', sprintf('GRFs_Fit_%02d.mat', 64)));
    screenSize = [1680 1050];
    fixPoint = screenSize/2;
    stimCenter = round(fixPoint+(52*[2 3]));
    for n=1:64
        RFs(n).centerposy = screenSize(2) - RFs(n).centerposy;
    end
    caccept = ones(1, length(RFs));
elseif strcmp(name, 'Isis')
    load(fullfile('/mnt/hpx/projects/MWNaturalPredict/Cem/Analysis/RFmaps/', sprintf('GRFs_Fit_v1.mat')))
    for n=1:32
        RFss(n).centerposx = RFs{n}.gauss2D.centerposx;
        RFss(n).centerposy = RFs{n}.gauss2D.centerposy;
        RFss(n).sigmaX = real(RFs{n}.gauss2D.sigmaX);
        RFss(n).sigmaY = double(RFs{n}.gauss2D.sigmaY);
        RFss(n).angle = RFs{n}.gauss2D.angle;
        RFss(n).label = sprintf('ch-%02d', n);
    end
    RFs = RFss;
    screenSize = [1680 1050];
    fixPoint = screenSize/2;
    stimCenter = [960+16 660-136];
    caccept = ~ismember(1:32, [9 13 14 16 22 23 24 26 31 32]);
elseif strcmp(name, 'Ares')
    load(fullfile('/mnt/hpx/projects/MWNaturalPredict/Cem/Analysis/RFmaps/', sprintf('GRFs_Fit_v1_ares.mat')))
    for n=1:32
        RFss(n).centerposx = RFs{n}.gauss2D.centerposx;
        RFss(n).centerposy = RFs{n}.gauss2D.centerposy;
        RFss(n).sigmaX = real(RFs{n}.gauss2D.sigmaX);
        RFss(n).sigmaY = double(RFs{n}.gauss2D.sigmaY);
        RFss(n).angle = RFs{n}.gauss2D.angle;
        RFss(n).label = sprintf('ch-%02d', n);
    end
    RFs = RFss;
    screenSize = [1680 1050];
    fixPoint = screenSize/2;
    stimCenter = [980 520];
    caccept = ismember(1:32, [1 14 17 23 27 28 29 30 31]);
end

% get RFs on Im
fx = [RFs.centerposx] - fixPoint(1);
fy = [RFs.centerposy] - fixPoint(2);

% get RF positions wrt stimcenter
stimShift = stimCenter - fixPoint;
sx = round(fx - stimShift(1)) + [300 ];
sy = round(fy - stimShift(2)) + [300 ];

statLib = {highStats, lowStats};
statTag = {'highStats', 'lowStats'};

newSave = fullfile(saveDir, name);
if ~exist(newSave, 'dir')
    mkdir(newSave);
end
for ch=1:length(RFs)
    if caccept(ch)
        cfg = []; n = 0;
        for st=1:length(statLib)
            thisStat = statLib{st};
            for ii=1:length(thisStat)
                n = n + 1;
                imDir = fullfile(setDir, thisStat(ii).originalName);
                cfg{n}.stimDir= stimDir;
                cfg{n}.stimName= thisStat(ii).name;
                cfg{n}.originalName = thisStat(ii).originalName;
                cfg{n}.statType = statTag{st};
                cfg{n}.ch = ch;
                cfg{n}.name = name;
                
                % important stuff
                cfg{n}.inputfile = imDir;
                cfg{n}.outputfile = newSave;
                cfg{n}.centerx = sx(ch);
                cfg{n}.centery = sy(ch);
                cfg{n}.cSize = 64;
                cfg{n}.sSize = 128;
                cfg{n}.ssimSize = 64;
                cfg{n}.make_movie = false;
                cfg{n}.make_image = false;
            end
        end
        outfiles = tdt_inpainting_AP(cfg, 'slurm');
    end
end

%% merge all channels
allOut = [];
for ii=1:length(outfiles)
    ESIload(outfiles{ii});
    allOut = [allOut out];
end
save(fullfile(saveDir, sprintf('allChOrient_%s.mat', name)), 'allOut');

