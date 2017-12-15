function data = merge_processed_files(files, savename, cleanupFlag)
% Merge single channel processed LFP/MUA files
addpath('/mnt/hpx/opt/ESIsoftware/matlab/')

if ~exist('files', 'var')
    if ~iscell(filename); return; end % exit if cancel
    files = cellfun(@(x) fullfile(pathname, x), filename, 'UniformOutput', false);
end

if ~exist('cleanupFlag', 'var'); cleanupFlag = true; end

fprintf('Load %s\n', files{1})
tok = strsplit(files{1}, '.');
ESIload(files{1},'-mat');
% firstFile = struct2cell(load(files{1}, '-mat'));
eval(sprintf('firstFile = %s;', tok{end}));
firstFile = firstFile{1};

nSamples = length(firstFile.data);
nChannels = length(files);

data = [];
data.trial = {zeros(nChannels, nSamples)};
data.trial{1}(1,:) = firstFile.data;
data.time = {(0:nSamples-1)/firstFile.fsample};
data.fsample = firstFile.fsample;
data.cfg = firstFile.cfg;
data.sampleinfo = [1 nSamples];

% channel labels
data.label = cell(nChannels,1);
if any(firstFile.header.channelNum) %not zero as in Atos
    data.label{1} = sprintf('%s_ch%04d', firstFile.datatype, firstFile.header.channelNum);
else
    [~,fname,~]=fileparts(files{1});
    tok = strsplit(fname, '_');
    data.label{1} = sprintf('%s_%s', firstFile.datatype, tok{end});
end

clear firstFile

for iFile = 2:nChannels
    fprintf('Load %s\n', files{iFile});
    [~,fname,~]=fileparts(files{iFile});
    tok = strsplit(fname, '_');
    chn = tok{end};
    chNum = str2num(chn(3:end));
    tok = strsplit(files{1}, '.');
    (ESIload(files{iFile}, '-mat'));
    eval(sprintf('currentFile = %s;', tok{end}));
    currentFile = currentFile{1};
    
    data.trial{1}(chNum,:) = currentFile.data;
    if any(currentFile.header.channelNum) %not zero as in Atos
        data.label{chNum} = sprintf('%s_ch%04d', currentFile.datatype, currentFile.header.channelNum);
    else
        [~,fname,~]=fileparts(files{iFile});
        tok = strsplit(fname, '_');
        data.label{chNum} = sprintf('%s_%s', currentFile.datatype, tok{end});
    end
end

if exist('savename', 'var') && ischar(savename)
    fprintf('Save to %s\n', savename);
    save(savename, 'data', '-v7.3')
end

if cleanupFlag
    fprintf('Clean up single channel files\n')
    for iFile = 1:length(files)
        delete(files{iFile})
    end
end

% clean up if no output arguments
if nargout == 0
    clear data
end