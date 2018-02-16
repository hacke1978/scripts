function tdt_analysis_AP(allCfg, calcLocation)
addpath /mnt/hpx/opt/ESIsoftware/matlab/;
addpath /mnt/hpx/opt/ESIsoftware/slurmfun/;

%LFP and MUA are slurmed
ths = [10 35 36 30 32 15 34 37 54 55 60 73:84];


switch calcLocation
    case 'slurm'
        license('inuse')
        slurmfun(@analyze_session, allCfg, 'stopOnError', false, 'partition', '16GBL', 'useUserPath', true, 'waitForToolboxes', {'statistics_toolbox', 'signal_toolbox'});
%         slurmfun(@analyze_session, allCfg, 'stopOnError', false, 'partition', '24GBL', 'useUserPath', true, 'waitForToolboxes', {'statistics_toolbox', 'signal_toolbox'});
%           slurmfun(@analyze_session, allCfg, 'stopOnError', false, 'partition', '48GBL', 'useUserPath', true, 'waitForToolboxes', {'statistics_toolbox', 'signal_toolbox'});
    case 'local'
        cellfun(@analyze_session, allCfg, 'UniformOutput', false);
end