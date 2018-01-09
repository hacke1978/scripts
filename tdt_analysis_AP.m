function tdt_analysis_AP(allCfg, calcLocation)
addpath /mnt/hpx/opt/ESIsoftware/matlab/;
addpath /mnt/hpx/opt/ESIsoftware/slurmfun/;
%LFP and MUA are slurmed

switch calcLocation
    case 'slurm'
%         originalDirectory = pwd();
%         cd(fullfile('/mnt/hpx/slurm/', getenv('USER'))); % 8GBS
        license('inuse')
        slurmfun(@analyze_session, allCfg, 'stopOnError', false, 'partition', '24GBL', 'useUserPath', true, 'waitForToolboxes', {'statistics_toolbox', 'signal_toolbox'});
%         cd(originalDirectory);
    case 'local'
        cellfun(@analyze_session, allCfg, 'UniformOutput', false);
end
