function data = merge_RF_files(cfg, savename, cleanupFlag)
% Merge single channel processed LFP/MUA files
addpath('/mnt/hpx/opt/ESIsoftware/matlab/')

% if ~exist('files', 'var')
%     if ~iscell(filename); return; end % exit if cancel
%     files = cellfun(@(x) fullfile(pathname, x), filename, 'UniformOutput', false);
% end

% if ~exist('cleanupFlag', 'var'); cleanupFlag = true; end
files = dir(fullfile(savename, 'ch*.mat'));

RFs = [];
for ii=1:length(files)
        fname = files(ii).name;
        tok = strsplit(fname, '_');
        ch = str2num(tok{1}(3:end));
        load(fullfile(savename, fname));
        RFs(ch).centerposx = RF.centerposx;
        RFs(ch).centerposy = cfg.screenSize(2)-RF.centerposy;
        RFs(ch).sigmaX = RF.sigmaX;
        RFs(ch).sigmaY = RF.sigmaY;
        RFs(ch).angle = RF.angle;
        RFs(ch).label_tdt = sprintf('ch-%02d', ch);
end

for ch=1:cfg.nChan
   if isempty(RFs(ch).centerposx)
        RFs(ch).centerposx = nan;
        RFs(ch).centerposy = nan;
        RFs(ch).sigmaX = nan;
        RFs(ch).sigmaY = nan;
        RFs(ch).angle = nan;
        RFs(ch).label_tdt = sprintf('ch-%02d', ch);
   end
end

% save 
tok = strsplit(savename, '/');
save(fullfile(savename, sprintf('%s.RFs', tok{end})), 'RFs')
% load(fullfile(savename, sprintf('%s.RFs', tok{end})), '-mat')
% load(fullfile('/mnt/hpx/projects/MWNaturalPredict/Hermes/rfmapping-bar', sprintf('%s.RF, tok{end}))), '-mat')
cfg = [];
cfg.orient = 0:45:360-45;
cfg.nChan = 128;
cfg.screenSize = [1680 1050];
cfg.fixPoint = cfg.screenSize/2;
% plot the rfs
fullScreen = ones(cfg.screenSize(2), cfg.screenSize(1))*128;
% h = figure; %set(h, 'visible', 'off');
% imagesc(fullScreen); colormap gray; hold on;
% plot(cfg.fixPoint(1), cfg.fixPoint(2), 'ro')
% RFlist = [32 10 16 14 23 41 37 35]
% for ch=1:length(RFs)
%     if ismember(ch, RFlist)
%    ellipsedraw(RFs(ch).sigmaX/2, RFs(ch).sigmaY/2, ...
%                 RFs(ch).centerposx, RFs(ch).centerposy, ...
%                 -RFs(ch).angle, 'k', [128 128], 0);hold on;
% %             text(RFs(ch).centerposx,RFs(ch).centerposy, data.label(ch), 'FontSize', 8, 'FontWeight', 'bold');
% text(RFs(ch).centerposx,RFs(ch).centerposy, {ch}, 'FontSize', 8, 'FontWeight', 'bold');
%     end
% end

allPoints = repmat([934 830] - [cfg.fixPoint], length(cfg.orient), 1);
cfg.orient = cfg.orient + 90;
(sum([allPoints(:, 1).*sin((cfg.orient)*pi/180)' -allPoints(:, 2).*cos((cfg.orient)*pi/180)'], 2)/1981*3300);

% % fprintf('Load %s\n', files{1})
% tok = strsplit(files{1}, '.');
% ESIload(files{1},'-mat');
% % firstFile = struct2cell(load(files{1}, '-mat'));
% eval(sprintf('firstFile = %s;', tok{end}));
% firstFile = firstFile{1};
% 
% nSamples = length(firstFile.data);
% nChannels = length(files);
% 
% data = [];
% data.trial = {zeros(nChannels, nSamples)};
% data.trial{1}(1,:) = firstFile.data;
% data.time = {(1:nSamples)/firstFile.fsample};
% data.fsample = firstFile.fsample;
% data.cfg = firstFile.cfg;
% data.sampleinfo = [1 nSamples];
% 
% % channel labels
% data.label = cell(nChannels,1);
% if any(firstFile.header.channelNum) %not zero as in Atos
%     data.label{1} = sprintf('%s-ch%03d', firstFile.datatype, firstFile.header.channelNum);
% else
%     [~,fname,~]=fileparts(files{1});
%     tok = strsplit(fname, '_');
%     chn = tok{end};
%     data.label{1} = sprintf('%s-%s', firstFile.datatype, chn(3:end));
% end
% 
% clear firstFile
% 
% for iFile = 2:nChannels
%     fprintf('Load %s\n', files{iFile});
%     [~,fname,~]=fileparts(files{iFile});
%     tok = strsplit(fname, '_');
%     chn = tok{end};
%     chNum = str2num(chn(3:end));
%     tok = strsplit(files{1}, '.');
%     (ESIload(files{iFile}, '-mat'));
%     eval(sprintf('currentFile = %s;', tok{end}));
%     currentFile = currentFile{1};
%     
%     data.trial{1}(chNum,:) = currentFile.data;
%     if any(currentFile.header.channelNum) %not zero as in Atos
%         data.label{chNum} = sprintf('%s-%03d', currentFile.datatype, currentFile.header.channelNum);
%     else
%         [~,fname,~]=fileparts(files{iFile});
%         tok = strsplit(fname, '_');
%         chn = tok{end};
%         data.label{chNum} = sprintf('%s-%s', currentFile.datatype, chn(3:end));
%     end
% end
% 
% if exist('savename', 'var') && ischar(savename)
%     fprintf('Save to %s\n', savename);
%     save(savename, 'data', '-v7.3')
% end
% 
% if cleanupFlag
%     fprintf('Clean up single channel files\n')
%     for iFile = 1:length(files)
%         delete(files{iFile})
%     end
% end
% 
% % clean up if no output arguments
% if nargout == 0
%     clear data
% end