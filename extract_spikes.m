function [outFile] = extract_spikes(cfg)
% Create LFP and MUA files from SEV files
addpath('/opt/fieldtrip/')
addpath('/mnt/hpx/opt/ESIsoftware/matlab/')
addpath('/opt/ESIsoftware/slurmfun/')
addpath('/mnt/hpx/slurm/uranc/TDT2MAT/TDT_MichaelStephan_Add2Path/')
file = cfg.filename;
% filename
if ~exist('file', 'var')
    %     [filename, pathname] = uigetfile(...
    %         {'*.sev'}, 'Pick a SEV file', '');
    if ~ischar(filename); return; end % exit if cancel
    file = fullfile(pathname, filename);
    clear filename pathname
end

[folder, basename, ext] = fileparts(file);

% configuration options
if ~exist('cfg', 'var'); cfg = []; end
if ~isfield(cfg, 'processSpike');    cfg.processSpike = true;    end
if ~isfield(cfg, 'targetFolder');  cfg.targetFolder = [];        end
if isempty(cfg.targetFolder);      cfg.targetFolder = folder;    end
if ~isfield(cfg, 'muaFreqBand');   cfg.muaFreqBand = [300 6000]; end
if ~isfield(cfg, 'pSigma');   cfg.pSigma = 3; end
if ~isfield(cfg, 'pISI');   cfg.pISI = 1.5; end
if ~isfield(cfg, 'save');          cfg.save = true;              end

spikeFile = '';

%% Read data
fprintf('Load %s ...', [basename ext])
[data, header] = read_data(file);
fprintf('\n')

if ~isfield(header,'Fs')
    header.Fs=2.441406250000000e+04; %fix for Atos, weird header.
end

header.channelNum = 0;
% header.channelNum=32; %fix for Atos, weird header.
%% MUA
muat = [];
if cfg.processSpike
    fprintf('Process MUA: ')
    
    muat.header = header;
    muat.cfg = cfg;
    muat.cfg.file = file;
    muat.datatype = 'muat';
    muat.sampleinfo = [1 length(data)];
    muat.data = double(data);
    muat.fsample = header.Fs;
    % construct butterworth band-pass filter
    [b,a] = butter(4, cfg.muaFreqBand/(header.Fs/2), 'bandpass');%butter(4, cfg.muaFreqBand/(header.Fs/2), 'bandpass');
    
    fprintf('filter & rectify')
    
    % filter data with zero-phase
    muat.data = filtfilt(b, a, muat.data);
    disp(cfg.pSigma)
    minThreshold =cfg.pSigma*median(abs(muat.data)/0.6745);
    %     minThreshold =3*std(muat.data./0.6745);
    minISI = round(cfg.pISI*1e-3*header.Fs);
    spikes = peakseek(muat.data, minISI, minThreshold);
    %     spikeTrain = zeros(size(mua.data));
    %     spikeTrain(spikes) = 1;
    muat.data = spikes/muat.fsample;
    
    muat = {muat};
    % save
    if cfg.save
        spikeFile = fullfile(cfg.targetFolder, [basename, '.muat']);
        fprintf(', save %s', spikeFile);
        try
            ESIsave(spikeFile, 'muat')
        catch me
            warning('Using -v7.3 flag for saving due to data length')
            ESIsave(spikeFile, 'muat')
        end
    end
    fprintf('\n')
end

outFile.spikeFile = spikeFile;
% clean up if no output arguments
if nargout == 0
    clear lfp mua data
end

