function tdt_figures_AP(allCfg, calcLocation)
addpath('/mnt/hpx/opt/ESIsoftware/matlab/')
addpath /mnt/hpx/opt/ESIsoftware/slurmfun/

%LFP and MUA are slurmed
switch calcLocation
    case 'slurm'
        originalDirectory = pwd();
        cd(fullfile('/mnt/hpx/slurm/', getenv('USER')));
        license('inuse')
        out = slurmfun(@plot_session, allCfg, 'partition', '8GBS', 'useUserPath', true, 'waitForToolboxes', {'statistics_toolbox'});
        cd(originalDirectory);
    case 'local'
        out = cellfun(@plot_session, allCfg, 'UniformOutput', false);
end

% plot group stuff
% corrCfg = [];
% corrCfg.outputfile = '/mnt/hpx/projects/MWNaturalPredict/Cem/Analysis/NatImFix';
% corrCfg.corrLayout = 'group';
% corrCfg.print = true;
% corrCfg.do_lfpPower2 = allCfg{1}.do_lfpPower2;
% corrCfg.isOri = false;
% plotCorrelations(corrCfg, out);
% 
% allOut = cat(1, out{:});
% poolFiles = cat(2, allOut.allFiles);
% catCfg = cat(1, allOut.allCfg);
% [baseName, ~] = fileparts(catCfg(1).outputfile);
% sesNames = cellfun(@(x) (x{end}), cellfun(@(x) strsplit(x, '/'), ...
%     {catCfg.outputfile}, 'UniformOutput', false), 'UniformOutput', false);
% poolCfg = allCfg{1};
% poolCfg.outputfile = fullfile(baseName, strcat(sesNames{:}));
% poolCfg.nConds = cellfun(@(x) length(x.allFiles), out);
% if ~exist(poolCfg.outputfile, 'dir'); mkdir(poolCfg.outputfile); end
% plot_all(poolCfg, poolFiles);
% plotPatchFromIm(poolCfg, poolFiles);