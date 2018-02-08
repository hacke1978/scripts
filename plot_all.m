function plot_all(allCfg, allFiles, varargin)
% Plot the analysis results
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Cem Uran
% cem.uran@esi-frankfurt.de
% 09/07/17
%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Choose layout according to setup
if strcmp(allCfg.name, 'Isis')
    layout = vertcat([NaN, NaN, 27, 1, NaN, NaN],[reshape(17:26, 5, 2) [28:32]' reshape(2:16, 5, 3)]);
elseif strcmp(allCfg.name, 'Hermes')
    layout = reshape(64:127, 8, 8);
elseif strcmp(allCfg.name, 'Ares')
    layout = vertcat([NaN, NaN, 27, 1, NaN, NaN],[reshape(17:26, 5, 2) [28:32]' reshape(2:16, 5, 3)]);
end

if nargin > 2
    overlayFiles = varargin{1};
end
% Start the figure
fileList = fieldnames(allFiles);

for fl = 1:length(fileList)
    thisFile = fileList{fl};
    if strcmp(thisFile, 'lfpPower')
        fprintf('plotting %s \n', thisFile)
        if allCfg.isOverlay
            plotLFPpower(allCfg, allFiles, thisFile, layout, overlayFiles)
        else
            plotLFPpower(allCfg, allFiles, thisFile, layout)
        end
        if (strcmp(thisFile, 'lfpPower') && allCfg.onIm && strcmp(allCfg.layout, 'channels'))
            plotOnIm(allCfg, allFiles, thisFile)
        end
    elseif ...
            ...%             strcmp(thisFile, 'muaxPPCSlepian') || ...
            ...%             strcmp(thisFile, 'muaxPPCHanning') || ...
            strcmp(thisFile, 'muaxCoherence') || ...
            strcmp(thisFile, 'timelockLFP') || ...
            strcmp(thisFile, 'timelockMUAX') || ...
            strcmp(thisFile, 'trialPSTH') || ...
            strcmp(thisFile, 'stSpec')
        fprintf('plotting %s \n', thisFile)
        if allCfg.isOverlay
            plot_single(allCfg, allFiles, thisFile, layout, overlayFiles);
        else
            plot_single(allCfg, allFiles, thisFile, layout);
        end
        if (strcmp(thisFile, 'stSpec') && ~allCfg.onIm)
            %             plotOnIm(allCfg, allFiles, thisFile)
        end
    else
        fprintf('NOT plotting %s \n', thisFile)
    end
end