function [out] = getPredForIm(filename)
% get the stats for im

tok  = strsplit(filename, '/');
if strcmp(tok{end-2}, '_stimset21')
    stimSetAll = '/mnt/v7k/home/uranc/workspace/VinckLab/Dataset/Processed/_6flickr/stimsetChosen/_stimset21';
    stripNo = strsplit(filename, '_');
    sesNo = str2num(stripNo{end-1}(1:2));
    imNo = str2num(stripNo{end-2}(end-1:end));
    stimsetDir = fullfile(stimSetAll, sprintf('stimset%02d', sesNo));
    load(fullfile(stimsetDir, 'chosenInd.mat'));
    load(fullfile(stimSetAll, 'allStats_CagliF.mat'));
    
    if sesNo < 5
        allInd = [chosenInd.highInd chosenInd.lowInd 0 0];
    else
        allInd = [chosenInd.highInd chosenInd.lowInd 0];
    end
    try
        randInd = chosenInd.randName;
    catch
        randInd = chosenInd.randInd;
    end
    
    if allInd(find(imNo==randInd)) == 0
        pred = 1;
        sname = 'gray';
    else
        pred = allStats(allInd(find(imNo==randInd))).ssimNaive(5);
        sname = allStats(allInd(find(imNo==randInd))).name;
    end
    out.pred = pred;
    out.sname = sname;
elseif strcmp(tok{end-2}, 'Andreea')
    stimSetAll = '/mnt/v7k/home/uranc/workspace/VinckLab/Dataset/Processed/_6flickr/stimsetChosen/Andreea';
    stripNo = strsplit(filename, '_');
    sesNo = str2num(stripNo{end-1}(1:2));
    imNo = str2num(stripNo{end-2}(end-1:end));
    if sesNo>17
        stimsetDir = fullfile(stimSetAll, sprintf('stimsetSelect%02d', sesNo));
    else
        stimsetDir = fullfile(stimSetAll, sprintf('stimsetNo%02d', sesNo));
    end
    load(fullfile(stimsetDir, 'chosenInd.mat'));
    load(fullfile(stimSetAll, 'allStats_CagliF.mat'));
    
    % Get Im
    orInd = cellfun(@(x) str2num(x(10:13)), vertcat(chosenInd.highIndSim, chosenInd.lowIndStim), 'UniformOutput', false);
    allInd = [horzcat(orInd{:})    0];
    try
        randInd = chosenInd.randName;
    catch
        randInd = chosenInd.randInd;
    end
    
    if allInd(find(imNo==randInd)) == 0
        pred = 1;
        sname = 'gray';
    else
        pred = allStats(allInd(find(imNo==randInd))).ssimNaive(5);
        sname = allStats(allInd(find(imNo==randInd))).name;
    end
    	out.pred = pred;
        out.sname = sname;
end