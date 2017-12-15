function [varargout] = getImListForSession(allCfg, dataDirSes)

if strcmp(allCfg.name, 'Isis')
    logPath = dir(fullfile(dataDirSes, '*_MatLog.mat'));
    load(fullfile(dataDirSes, logPath.name));
    stimList = vertcat({Config.Protocol.ImageList.name});
    
    stimsetAll = '/mnt/v7k/home/uranc/workspace/VinckLab/Dataset/Processed/_6flickr/stimsetChosen/Andreea';
    sesNo = []; stimDir = '';
    for is = 1:length(stimList)
        stimParam = strsplit(stimList{is}, '_');
        sesNo = str2num(stimParam{2});
        
        if sesNo>17
            stimDir(is, :) = fullfile(stimsetAll, sprintf('stimsetSelect%02d', sesNo));
        elseif strcmp('/mnt/hpx/projects/MWNaturalPredict/Isis/isisc472a01', dataDirSes)
            stimDir(is, :) = fullfile(stimsetAll, sprintf('stimsetTuning%02d', sesNo));
        else
            stimDir(is, :) = fullfile(stimsetAll, sprintf('stimsetNo%02d', sesNo));
        end
    end
    imList = strcat(stimDir, '/', stimList');
    varargout{1} = imList;
elseif strcmp(allCfg.name, 'Hermes')
    imNames = dir(fullfile(dataDirSes, '*.png'));
    imList = strcat(dataDirSes, '/', {imNames.name}');
    varargout{1} = imList;
elseif strcmp(allCfg.name, 'Ares')
    logPath = dir(fullfile(dataDirSes, '*_MatLog.mat'));
    load(fullfile(dataDirSes, logPath.name));
    stimsetAll = '/mnt/v7k/home/uranc/workspace/VinckLab/Dataset/Processed/_6flickr/stimsetChosen';
    if strcmp(allCfg.type, 'NatImFix')
        stimList = vertcat({Config.Protocol.ImageList.name});
        sesNo = []; stimDir = '';
        for is = 1:length(stimList)
            stimParam = strsplit(stimList{is}, '_');
            sesNo = str2num(stimParam{2});
            if sesNo > 67
                stimDir(is, :) = fullfile(stimsetAll, sprintf('Ares_stimsetSEQColorSizeTuning%02d', sesNo));
            else
                stimDir(is, :) = fullfile(stimsetAll, sprintf('Ares_stimsetSEQNatIm%02d', sesNo));
            end
        end
        imList = strcat(stimDir, '/', stimList');
        %         imNames = dir(fullfile(dataDirSes, '*.png'));
        %         imList = strcat(dataDirSes, '/', {imNames.name}');
        varargout{1} = imList;
    elseif strcmp(allCfg.type, 'NatImSEQ')
        stim = vertcat({Config.Protocol.Condition.CondFolder})';
        stim = cellfun(@(x) strsplit(x, '\'), stim, 'UniformOutput', false);
        stimDir = vertcat(stim{:});
        stimCond = stimDir(:, end);
        stimDir = stimDir{1, end-1};
        %         stimsetAll = '/mnt/v7k/home/uranc/workspace/VinckLab/Dataset/Processed/_6flickr/stimsetChosen';
        stimList = cellfun(@(x) dir(fullfile(stimsetAll, stimDir, x, '*.png')), stimCond, 'UniformOutput', false);
        stimList = cat(2, stimList{:});
        stimList = vertcat({stimList(1, :).name});
        %     stimList = vertcat({stimList(:, 1).name});
        
        % sort condition numbers 1 to end
        scond = cellfun(@(x) str2num(sprintf('%02d', str2num(x))), stimCond, 'UniformOutput', false);
        [~, sortInd] = sort([scond{:}]);
        imList = strcat(stimsetAll, '/', stimDir, '/', stimCond, '/', stimList');
        varargout{1} = imList;
        varargout{2} = sortInd;
    end
end