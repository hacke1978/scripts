function tdt_analysis_AP(allCfg, calcLocation)
addpath /mnt/hpx/opt/ESIsoftware/matlab/;
addpath /mnt/hpx/opt/ESIsoftware/slurmfun/;

%LFP and MUA are slurmed
switch calcLocation
    case 'slurm'
        license('inuse')
        slurmfun(@analyze_session, allCfg, 'stopOnError', false, 'partition', '16GBXL', 'useUserPath', true, 'waitForToolboxes', {'statistics_toolbox', 'signal_toolbox'});
%           slurmfun(@analyze_session, allCfg, 'stopOnError', false, 'partition', '48GBL', 'useUserPath', true, 'waitForToolboxes', {'statistics_toolbox', 'signal_toolbox'});
    case 'local'
        cellfun(@analyze_session, allCfg, 'UniformOutput', false);
end