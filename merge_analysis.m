function data = merge_analysis(files, cleanupFlag)
% Merge single channel processed LFP/MUA files
addpath('/mnt/hpx/opt/ESIsoftware/matlab/')

if ~exist('cleanupFlag', 'var'); cleanupFlag = true; end


outFiles = cat(2, files{:});

for fl=fieldnames(outFiles)'
    thisList = cat(2, outFiles.(fl{1}));
    for tr=1:size(thisList, 1)
        fprintf('Load %s\n', thisList{tr, 1});
        tok = strsplit(thisList{tr, 1}, '/');
        savename = [thisList{tr, 1}(1:strfind(thisList{tr, 1}, 'ch')-2) '.mat'];
        for ch=1:size(thisList, 2)
            if ch == 1
                ESIload(thisList{tr, ch},'-mat');
                if strcmp(fl{:}, 'staFiles')
                    all_stAv = stAv;
                elseif strcmp(fl{:}, 'sfcFiles')
                    all_stSpec = stSpec;
                elseif strcmp(fl{:}, 'sfcFilesVar')
                    all_stSpecVar = stSpecPerTrial;
                end
            else
                ESIload(thisList{tr, ch},'-mat');
                if strcmp(fl{:}, 'staFiles')
                    all_stAv = vertcat(all_stAv, stAv);
                elseif strcmp(fl{:}, 'sfcFiles')
                    all_stSpec = vertcat(all_stSpec, stSpec);
                elseif strcmp(fl{:}, 'sfcFilesVar')
                    all_stSpecVar = vertcat(all_stSpecVar, stSpecPerTrial);
                end
            end
        end
        % save and clean up
        if strcmp(fl{:}, 'staFiles')
            stAv = all_stAv;
            ESIsave(savename, 'stAv')
            clear all_stAv stAv
        elseif strcmp(fl{:}, 'sfcFiles')
            stSpec = all_stSpec;
            ESIsave(savename, 'stSpec')
            clear all_stSpec stSpec
        elseif strcmp(fl{:}, 'sfcFilesVar')
            stSpecPerTrial = all_stSpecVar;
            ESIsave(savename, 'stSpecPerTrial')
            clear all_stSpecVar stSpecPerTrial
        end
    end
    fprintf('Clean up single channel files\n')
    for tFile = 1:size(thisList, 1)
        for cFile = 1:size(thisList, 2)
            delete(thisList{tFile, cFile})
        end
    end
end

% clean up if no output arguments
if nargout == 0
    clear data
end
