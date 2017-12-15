function data = merge_spike_files(files, savename, cleanupFlag)
% Merge single channel processed LFP/MUA files
addpath('/mnt/hpx/opt/ESIsoftware/matlab/')
if ~exist('files', 'var')
    [filename, pathname] = uigetfile(...
        {'*.mua;*.lfp', 'MUA/LFP files (*.mua, *.lfp)'; ...
        '*.mua', 'MUA files (*.mua)'; ...
        '*.lfp', 'LFP files (*.lfp)'}, ...
        'Pick a MUA/LFP file', '', ...
        'MultiSelect', 'on');
    if ~iscell(filename); return; end % exit if cancel
    files = cellfun(@(x) fullfile(pathname, x), filename, 'UniformOutput', false);
end

if ~exist('cleanupFlag', 'var'); cleanupFlag = true; end

fprintf('Load %s\n', files{1})
% firstFile = struct2cell(load(files{1}, '-mat'));
tok = strsplit(files{1}, '.');
ESIload(files{1},'-mat');
eval(sprintf('firstFile = %s;', tok{end}));
firstFile = firstFile{1};

nSamples = length(firstFile.data);
nChannels = length(files);

data = [];
data.trial = {zeros(nChannels, nSamples)};
data.trial{1} = firstFile.data;
data.fsample = firstFile.fsample;
data.sampleinfo = firstFile.sampleinfo;
data.cfg = firstFile.cfg;
data.hdr = firstFile.header;
% channel labels
data.label = cell(nChannels,1);
if any(firstFile.header.channelNum) %not zero as in Atos
    data.label{1} = sprintf('%s-%04d', firstFile.datatype, firstFile.header.channelNum);
else
    [~,fname,~]=fileparts(files{1});
    tok = strsplit(fname, '_');
    chn = tok{end};
    data.label{1} = sprintf('%s-%s', firstFile.datatype, chn(3:end));
end

clear firstFile
for iFile = 2:nChannels
    fprintf('Load %s\n', files{iFile});
%     currentFile = struct2cell(load(files{iFile}, '-mat'));
    [~,fname,~]=fileparts(files{iFile});
    tok = strsplit(fname, '_');
    chn = tok{end};
    chNum = str2num(chn(3:end));
    tok = strsplit(files{1}, '.');
    ESIload(files{iFile}, '-mat');
    eval(sprintf('currentFile = %s;', tok{end}));
    currentFile = currentFile{1};
    
    data.trial{chNum} = currentFile.data;
    if any(currentFile.header.channelNum) %not zero as in Atos
        data.label{chNum} = sprintf('%s-ch%04d', currentFile.datatype, currentFile.header.channelNum);
    else
        [~,fname,~]=fileparts(files{iFile});
        tok = strsplit(fname, '_');
        chn = tok{end};
        data.label{chNum} = sprintf('%s-%s', currentFile.datatype, chn(3:end));
    end
end

if exist('savename', 'var') && ischar(savename)
    fprintf('Save to %s\n', savename);
    try
        save(savename, 'data', '-v7.3')
    catch me
        warning('Using -v7.3 flag for saving due to data length')
        save(savename, 'data', '-v7.3')
    end
    
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
