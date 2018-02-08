function tdt_preprocessing_AP(cfg)
addpath('/opt/ESIsoftware/slurmfun/')
if ~exist('cfg', 'var'); cfg = []; end

if ~isfield(cfg, 'filename')
    [filename, pathname] = uigetfile(...
        {[ path '/*xWav*.sev']}, 'Pick raw SEV files', '', 'MultiSelect', 'on');
    filename = {filename};
    if ~iscell(filename); return; end % exit if cancel
    filename = cellfun(@(x) fullfile(pathname, x), filename, ...
        'UniformOutput', false);
else
    filename=cfg.filename;
end

if ~isfield(cfg, 'calcLocation'); cfg.calcLocation = 'slurm'; end % slurm | local
if ~isfield(cfg, 'processLfp');    cfg.processLfp = true;        end
if ~isfield(cfg, 'processMua');    cfg.processMua = true;        end
if ~isfield(cfg, 'processEms');    cfg.processEms = true;        end
if ~isfield(cfg, 'processDio');    cfg.processDio = true;        end
if ~isfield(cfg, 'decimation');    cfg.decimation = [8 3];      end %8 3 atos to get down to hermes/klecks. 4 4 to get down to previous atos sr.
if ~isfield(cfg, 'targetFolder');  cfg.targetFolder = [];        end
if isempty(cfg.targetFolder);      cfg.targetFolder = pathname;  end
if ~isfield(cfg, 'muaFreqBand');   cfg.muaFreqBand  = [300 12000]; end
if ~isfield(cfg, 'save');          cfg.save = true;              end

allCfg = repmat({cfg}, 1, length(filename));
for ii=1:length(filename)
    allCfg{ii}.filename = filename{ii};
end

% get events and transorm into BR style
% LFP and MUA are slurmed

switch cfg.calcLocation
    case 'slurm'
%         originalDirectory = pwd();
%         cd(fullfile('/mnt/hpx/slurm/', getenv('USER')));
        license('inuse')
        out = slurmfun(@make_lfp_mua, allCfg, 'partition', '16GBL',  'waitForToolboxes', {'signal_toolbox'}, 'stopOnError', false);
%         cd(originalDirectory);
    case 'local'
        out = cellfun(@make_lfp_mua, allCfg, 'UniformOutput', false);
end

outFiles = cat(2, out{:});
lfpFiles = cat(2, {outFiles.lfpFile});
muaFiles = cat(2, {outFiles.muaFile});
rawMuaFiles = cat(2, {outFiles.rawMuaFile});
lfpFile = [lfpFiles{1}(1:strfind(lfpFiles{1}, 'xWav')+3) '.lfp'];
merge_processed_files(lfpFiles, lfpFile, true)
muaFile = [muaFiles{1}(1:strfind(muaFiles{1}, 'xWav')+3) '.muax'];
merge_processed_files(muaFiles, muaFile, true)
rawMuaFile = [rawMuaFiles{1}(1:strfind(rawMuaFiles{1}, 'xWav')+3) '.mua'];
merge_processed_files(rawMuaFiles, rawMuaFile, true)
