function filedir = tdt_inpainting_AP(allCfg, calcLocation)
addpath('/mnt/hpx/opt/ESIsoftware/matlab/')
addpath /mnt/hpx/opt/ESIsoftware/slurmfun/
%LFP and MUA are slurmed

switch calcLocation
    case 'slurm'
        originalDirectory = pwd();
        cd(fullfile('/mnt/hpx/slurm/', getenv('USER')));
        license('inuse')
%         out = slurmfun(@getInpainting, allCfg, 'partition', '8GBS', 'useUserPath', true, 'waitForToolboxes', {'image_toolbox'});
        out = slurmfun(@getOrientation, allCfg, 'partition', '8GBS', 'useUserPath', true, 'waitForToolboxes', {'image_toolbox'});
        cd(originalDirectory);
    case 'local'
%         out = cellfun(@getInpainting, allCfg, 'UniformOutput', false);
        out = cellfun(@getOrientation, allCfg, 'UniformOutput', false);
end
ch = allCfg{1}.ch;
% st = allCfg{1}.statType;
savename = allCfg{1}.outputfile;
filedir = fullfile(savename, sprintf('ch%02d_orient.mat', ch));
ESIsave(filedir, 'out');
end