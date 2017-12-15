function tdt_extractspikes_AP(cfg)
addpath('/opt/ESIsoftware/slurmfun/')

if ~exist('cfg', 'var'); cfg = []; end

if ~isfield(cfg, 'filename')
    filename = {filename};
    if ~iscell(filename); return; end % exit if cancel
    filename = cellfun(@(x) fullfile(pathname, x), filename, ...
                        'UniformOutput', false);
else
    filename=cfg.filename;
end

if ~isfield(cfg, 'calcLocation'); cfg.calcLocation = 'slurm'; end % slurm | local
if ~isfield(cfg, 'processDio');    cfg.processSpike = true;        end
if ~isfield(cfg, 'targetFolder');  cfg.targetFolder = pathname;  end
if ~isfield(cfg, 'muaFreqBand');   cfg.muaFreqBand  = [300 12000]; end
if ~isfield(cfg, 'pSigma');   cfg.pSigma = 3; end
if ~isfield(cfg, 'pISI');   cfg.pISI = 1.5; end
if ~isfield(cfg, 'save');          cfg.save = true;              end

allCfg = repmat({cfg}, 1, length(filename));
for ii=1:length(filename)
    allCfg{ii}.filename = filename{ii};
end

%LFP and MUA are slurmed

switch cfg.calcLocation
    case 'slurm'
        originalDirectory = pwd();
        cd(fullfile('/mnt/hpx/slurm/', getenv('USER')));
        out = slurmfun(@extract_spikes, allCfg, 'partition', '8GBS', 'useUserPath', true, 'waitForToolboxes', {'signal_toolbox'});
        cd(originalDirectory);
    case 'local'
        out = cellfun(@extract_spikes, allCfg, 'UniformOutput', false);
end

outFiles = cat(2, out{:});
spikeFiles = cat(2, {outFiles.spikeFile});
spikeFile = [spikeFiles{1}(1:strfind(spikeFiles{1}, 'xWav')+3) '.muat'];
merge_spike_files(spikeFiles, spikeFile, true)

