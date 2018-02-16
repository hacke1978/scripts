function [outFile] = make_lfp_mua(cfg)
addpath('/opt/ESIsoftware/matlab/')
% Create LFP and MUA files from SEV files
file = cfg.filename;
[folder, basename, ext] = fileparts(file);

% configuration options
if ~exist('cfg', 'var'); cfg = []; end
if ~isfield(cfg, 'processLfp');    cfg.processLfp = true;        end
if ~isfield(cfg, 'processMua');    cfg.processMua = true;        end
if ~isfield(cfg, 'decimation');    cfg.decimation = [4 4];       end
if ~isfield(cfg, 'targetFolder');  cfg.targetFolder = [];        end
if isempty(cfg.targetFolder);      cfg.targetFolder = folder;    end
if ~isfield(cfg, 'muaFreqBand');   cfg.muaFreqBand = [300 6000]; end
if ~isfield(cfg, 'save');          cfg.save = true;              end

lfpFile = '';
muaFile = '';

%% Read data
fprintf('Load %s ...', [basename ext])
[data, header] = read_data(file);
fprintf('\n')

if ~isfield(header,'Fs')
    header.Fs=2.441406250000000e+04; %fix for Atos, weird header.
end

header.channelNum = 0;
% LFP
lfp = [];
if cfg.processLfp
    fprintf('Process LFP: ')
    
    lfp.hdr.TimeStampPerSample = prod(cfg.decimation);
    lfp.hdr.Fs = header.Fs;
    lfp.hdr.FirstTimeStamp = 0;
    lfp.header = header;
    lfp.cfg = cfg;
    lfp.cfg.file = file;
    lfp.datatype = 'lfp';
    lfp.data = double(data);
    fprintf('decimate')
    lfp.fsample = header.Fs/prod(cfg.decimation);
    for iDecimation = 1:length(cfg.decimation)
        lfp.data = decimate(lfp.data, cfg.decimation(iDecimation), 'FIR');
    end
    lfp.date = datestr(now());
    
    lfp = {lfp};
    % ESIsave
    if cfg.save
        lfpFile = fullfile(cfg.targetFolder, [basename, '.lfp']);
        fprintf(', save %s', lfpFile)
        try
            ESIsave(lfpFile, 'lfp')
        catch me
            %             warning('Using -v7.3 flag for saving due to data length')
            ESIsave(lfpFile, 'lfp')
        end
    end
    fprintf('\n')
end

%% MUA
muax = [];
if cfg.processLfp
    fprintf('Process MUA: ')
    
    muax.header = header;
    muax.cfg = cfg;
    muax.cfg.file = file;
    muax.datatype = 'muax';
    
    muax.data = double(data);
    
    % construct butterworth band-pass filter
    [b,a] = butter(4, cfg.muaFreqBand/(header.Fs/2), 'bandpass');%butter(4, cfg.muaFreqBand/(header.Fs/2), 'bandpass');
    
    fprintf('filter & rectify')
    
    % filter data with zero-phase
    muax.data = filtfilt(b, a, muax.data);
    
    % get raw mua here
    mua = muax;

    % rectify
    muax.data  = abs(muax.data);
    
    % decimate
    fprintf(', decimate')
    for iDecimation = 1:length(cfg.decimation)
        muax.data = decimate(muax.data, cfg.decimation(iDecimation), 'FIR');
        mua.data = decimate(mua.data, cfg.decimation(iDecimation), 'FIR');
    end
    muax.date = datestr(now());
    muax.fsample = header.Fs/prod(cfg.decimation);
    muax = {muax};
    
    mua.date = datestr(now());
    mua.fsample = header.Fs/prod(cfg.decimation);
    mua = {mua};
    % ESIsave
    if cfg.save
        muaFile = fullfile(cfg.targetFolder, [basename, '.muax']);
        rawMuaFile = fullfile(cfg.targetFolder, [basename, '.mua']);
        fprintf(', save %s', rawMuaFile)
        try
            ESIsave(muaFile, 'muax')
            ESIsave(rawMuaFile, 'mua')
        catch me
            warning('Using -v7.3 flag for saving due to data length')
            ESIsave(muaFile, 'muax')
            ESIsave(rawMuaFile, 'mua')
        end
    end
    fprintf('\n')
end
outFile.lfpFile = lfpFile;
outFile.muaFile = muaFile;
outFile.rawMuaFile = rawMuaFile;
% clean up if no output arguments
if nargout == 0
    clear lfp mua data
end

