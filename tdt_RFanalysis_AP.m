function tdt_RFanalysis_AP(allCfg, calcLocation)
addpath('/mnt/hpx/opt/ESIsoftware/matlab/')
addpath /mnt/hpx/opt/ESIsoftware/slurmfun/
%LFP and MUA are slurmed

switch calcLocation
    case 'slurm'
        originalDirectory = pwd();
        cd(fullfile('/mnt/hpx/slurm/', getenv('USER')));
        license('inuse')
        slurmfun(@analyze_RFmaps, allCfg, 'partition', '8GBS', 'useUserPath', true, 'waitForToolboxes', {'statistics_toolbox', 'signal_toolbox'});
        cd(originalDirectory);
    case 'local'
        cellfun(@analyze_RFmaps, allCfg, 'UniformOutput', false)
 end