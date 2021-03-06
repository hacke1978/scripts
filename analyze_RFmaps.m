function analyze_RFmaps(allCfg)
addpath('/mnt/hpx/slurm/uranc/fieldtrip/')
% Requires on Fiorani, Mario, et al.
%"Automatic mapping of visual cortex receptive fields: a fast and precise algorithm."
% Journal of neuroscience methods 221 (2014): 112-126.
% Author: Cem Uran (cem.uran@esi-frankfurt.de)

% Parameters
filename = allCfg.inputfile;
screenSize = [1680 1050];
fixPoint = screenSize/2;

if strcmp(allCfg.name, 'Hermes')
    % %     muaName = dir(fullfile(filename, sprintf('*.muax')));
    % %     load(fullfile(filename, muaName.name ), '-mat');
    % %     label = data.label;
    % %     label_tdt = data.label_tdt;
    % %     rfName = dir(fullfile(filename, sprintf('*.RF')));
    % %     load(fullfile(filename, rfName.name ), '-mat');
    % %     orient = 0:45:360-45;
    % %     caccept = strncmp(label, 'V1', 2);
    % %     RFs = RF(caccept);
    % %     fit_gaussian_psth(allCfg, RFs, screenSize, label(caccept), label_tdt(caccept))
    calcLocation = 'slurm';
    chChosen = 1:128;
    
    cfg = cell(length(chChosen), 1);
    for ii=1:length(chChosen)
        ch = chChosen(ii);
        cfg{ii}.filename = filename;
        cfg{ii}.savename = fullfile(allCfg.outputfile, sprintf('ch%02d', ch));
        cfg{ii}.tag = '';
        cfg{ii}.orient = 0:45:359;
        cfg{ii}.caccept = ch;
        cfg{ii}.binSize = 0.005;
        cfg{ii}.pxd = round(degDistances(1, allCfg.name));
    end
    
else
    calcLocation = 'slurm';
    chChosen = 1:32;
    
    cfg = cell(length(chChosen), 1);
    for ii=1:length(chChosen)
        ch = chChosen(ii);
        cfg{ii}.tag = '';
        cfg{ii}.filename = filename;
        cfg{ii}.savename = fullfile(allCfg.outputfile, sprintf('ch%02d', ch));
        cfg{ii}.orient = 0:22.5:360-22.5;
        cfg{ii}.caccept = ch;
        cfg{ii}.binSize = 0.010;
        cfg{ii}.pxd = round(degDistances(1, allCfg.name));
    end
end

% choose where to compute
switch calcLocation
    case 'slurm'
        license('inuse')
        slurmfun(@getRFfromMUA, cfg, 'partition', '8GBS', 'useUserPath', true, 'stopOnError', false, 'waitForToolboxes', {'statistics_toolbox', 'curve_fitting_toolbox', 'image_toolbox'});
    case 'local'
        cellfun(@getRFfromMUA, cfg, 'UniformOutput', false)
end

cfg = [];
cfg.orient = 0:45:360-45;
cfg.savename = allCfg.outputfile;
cfg.nChan = 128;
cfg.screenSize = [1680 1050];
cfg.fixPoint = cfg.screenSize/2;
merge_RF_files(cfg, allCfg.outputfile, false)

end


%%
