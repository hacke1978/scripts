function out = plot_session(allCfg)
% Plot the analysis results
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Cem Uran
% cem.uran@esi-frankfurt.de
% 09/07/17
%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Params
filename = allCfg.inputfile;
savefile = allCfg.outputfile;

% Load all data
allFiles = load_session(allCfg, filename);

%% RF Dir
if strcmp(allCfg.name, 'Hermes')
    % load/reshape RFs manually for now
    load(fullfile('/mnt/hpx/projects/MWNaturalPredict/Cem/Analysis/RFmaps/', sprintf('GRFs_Fit_%02d.mat', 64)));
    screenSize = [1680 1050];
    fixPoint = screenSize/2;
    caccept = ones(64, 1); caccept(7) = 0; caccept(33) = 0;
%     caccept = strncmp([allFiles(1).RFs.label], 'V1', 10);
    if strcmp(allCfg.type, 'grating-ori')
        stimCenter = round(fixPoint+(51*[2 5]));
    elseif strcmp(filename, 'hermes_20171201_fixation-naturalim_90')
        stimCenter = round(fixPoint+(degDistances(1, allCfg.name)*[2 6]));
    else
        stimCenter = round(fixPoint+(degDistances(1, allCfg.name)*[2 3]));
    end
    for ii=1:64
        RFs(ii).centerposy = screenSize(2) - RFs(ii).centerposy;
    end
    
    % get orientation for the RF
    if allCfg.isOri
        readOri = xlsread('/mnt/hpx/projects/MWNaturalPredict/hermesRecordingLogAa.xlsx', 'Channels');
        ori = readOri(:, 3); ch = readOri(:, 1);
        cmap = cellfun(@(x) find(strncmp(cellfun(@num2str, num2cell(ch), 'UniformOutput', false), ...
            x(4:end), 5)), vertcat(RFs.label), 'UniformOutput', false);
        ori = num2cell(ori([cmap{:}]));
        [RFs(:).ori] = deal(ori{:});
    end
    
    % get predictibility measures
    if allCfg.isPred
        [statOut] = (cellfun(@(x) getPredForIm(x), imList, 'UniformOutput', false));
        statOut = cat(2, statOut{:});
        allPred = vertcat({statOut.pred}); allSname = vertcat({statOut.sname});
        [allFiles.pred] = deal(allPred{:}); [allFiles.sname] = deal(allSname{:});
    end
    
    % get images
    tok = strsplit(filename, '_');
    dataDir = '/mnt/hpx/projects/MWNaturalPredict/Hermes/Stimuli';
    if strcmp(allCfg.type, 'grating-ori')
        fname = fullfile(dataDir, sprintf('grating-ori'))
    else
        if isempty(str2num(tok{end}))
            sesNo = str2num(tok{end}(1:end-1));
        else
            sesNo = str2num(tok{end});
        end
        fname = fullfile(dataDir, sprintf('stimset%02d', sesNo))
%         load('/mnt/v7k/projects/MWNaturalPredict/StimuliCondInfo/hermes_20171113_fixation-color-masksize_78_condInfo.mat');
%         imList = fullfile(fname, condInfo.imName);
        imList = getImListForSession(allCfg, fname);
        [allFiles.imName] = deal(imList{:});
    end
    allFiles(1).RFs = RFs;
    allFiles(1).stimCenter = stimCenter;
    allFiles(1).caccept = caccept;
    
elseif strcmp(allCfg.name, 'Isis')
    load(fullfile('/mnt/hpx/projects/MWNaturalPredict/Cem/Analysis/RFmaps/', sprintf('GRFs_Fit_v1.mat')))
    for ii=1:32
        RFss(ii).centerposx = RFs{ii}.gauss2D.centerposx;
        RFss(ii).centerposy = RFs{ii}.gauss2D.centerposy;
        RFss(ii).sigmaX = real(RFs{ii}.gauss2D.sigmaX);
        RFss(ii).sigmaY = double(RFs{ii}.gauss2D.sigmaY);
        RFss(ii).angle = RFs{ii}.gauss2D.angle;
        RFss(ii).label = sprintf('ch-%02d', ii);
    end
    RFs = RFss;
    stimCenter = [960+16 660-136];
    caccept = ~ismember(1:32, [9 13 14 16 22 23 24 26 31 32]);
    
    % get images
    tok = strsplit(filename, '/');
    dataDir = '/mnt/hpx/projects/MWNaturalPredict/Isis/';
    fname = fullfile(dataDir, tok{end})
    imList = getImListForSession(allCfg, fname);
    
    % load all
    [allFiles.imName] = deal(imList{:});
    allFiles(1).RFs = RFs;
    allFiles(1).stimCenter = stimCenter;
    allFiles(1).caccept = caccept;
    
    if allCfg.isPred
        [statOut] = (cellfun(@(x) getPredForIm(x), imList, 'UniformOutput', false));
        statOut = cat(2, statOut{:});
        allPred = vertcat({statOut.pred}); allSname = vertcat({statOut.sname});
        [allFiles.pred] = deal(allPred{:}); [allFiles.sname] = deal(allSname{:});
    end
elseif strcmp(allCfg.name, 'Ares')
    load(fullfile('/mnt/hpx/projects/MWNaturalPredict/Cem/Analysis/RFmaps/', sprintf('GRFs_Fit_v1_ares.mat')))
    for ii=1:32
        RFss(ii).centerposx = RFs{ii}.gauss2D.centerposx;
        RFss(ii).centerposy = RFs{ii}.gauss2D.centerposy;
        RFss(ii).sigmaX = real(RFs{ii}.gauss2D.sigmaX);
        RFss(ii).sigmaY = double(RFs{ii}.gauss2D.sigmaY);
        RFss(ii).angle = RFs{ii}.gauss2D.angle;
        RFss(ii).label = sprintf('ch-%02d', ii);
    end
    RFs = RFss;
    stimCenter = [980 520];
    caccept = ismember(1:32, [14 26 27 29 30 31 32]);
    % images
    tok = strsplit(filename, '/');
    dataDir = fullfile('/mnt/hpx/projects/MWNaturalPredict/', allCfg.name, allCfg.type);
    fname = fullfile(dataDir, tok{end})
    if strcmp(allCfg.type, 'NatImFix')
        imList = getImListForSession(allCfg, fname);
    elseif strcmp(allCfg.type, 'NatImSEQ')
        [imList, sortInd] = getImListForSession(allCfg, fname);
        allFiles(1).sortInd = sortInd;
    end
    [allFiles.imName] = deal(imList{:});
    allFiles(1).RFs = RFs;
    allFiles(1).stimCenter = stimCenter;
    allFiles(1).caccept = caccept;
    
    if allCfg.isPred
        [statOut] = (cellfun(@(x) getPredForIm(x), imList, 'UniformOutput', false));
        statOut = cat(2, statOut{:});
        allPred = vertcat({statOut.pred}); allSname = vertcat({statOut.sname});
        [allFiles.pred] = deal(allPred{:}); [allFiles.sname] = deal(allSname{:});
    end
end

%% Plot things
plot_all(allCfg, allFiles);
plotPatchFromIm(allCfg, allFiles);
out = [];
out.allFiles = allFiles;
out.allCfg = allCfg;
% if allCfg.isCorr
%     rates = getRatesAndStats(allCfg, allFiles);
%     plotCorrelations(allCfg, rates);
% end