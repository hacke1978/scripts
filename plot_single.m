function plot_single(allCfg, allFiles, thisFile, layout)
% plot stuff on the image
% should i deally will be plotting automatically
savefile = allCfg.outputfile;
allThis = squeeze(cat(3, allFiles.(thisFile)));
caccept = allFiles(1).caccept;
shaded = 0;

if strcmp(thisFile, 'muaxCoherence')
    tlabel = allThis(1).labelcmb;
    tx = allThis(1).freq;
    data = cat(3, allThis.cohspctrm);
    txlim = [min(tx) max(tx)];
elseif strcmp(thisFile, 'muaxPPCSlepian')
    tlabel = allThis(1).labelcmb;
    tx = allThis(1).freq;
    data = cat(3, allThis.ppcspctrm);
    txlim = [min(tx) max(tx)];
elseif strcmp(thisFile, 'muaxPPCHanning')
    tlabel = allThis(1).labelcmb;
    tx = allThis(1).freq;
    data = cat(3, allThis.ppcspectrum);
    txlim = [min(tx) max(tx)];
elseif strcmp(thisFile, 'stSpec')
    tlabel = allThis(1).labelcmb(:, 2);
    if strcmp(allCfg.name, 'Isis')
        tlabel = cellfun(@ (x) (x(7:end)), tlabel, 'UniformOutput', false);
    elseif strcmp(allCfg.name, 'Ares')
        tlabel = cellfun(@ (x) (x(7:end)), tlabel, 'UniformOutput', false);
    elseif strcmp(allCfg.name, 'Hermes')
        tlabel = cellfun(@ (x) (x(4:end)), tlabel, 'UniformOutput', false);
        tlabel(find(strncmp(tlabel, 'X', 1))) = {'127'};
    end
    tx = allThis(1).freq;
    if allCfg.withErrorBars
        stsVar = zeros([size(allFiles(1).stSpecPerTrial(1).ppc1) size(allFiles, 2) size(allFiles(1).stSpecPerTrial, 1)]);
        for ii=1:size(allFiles, 2)
            for ch=1:size(allFiles(1).stSpecPerTrial, 1)
                data_var = cat(3, allFiles(ii).stSpecPerTrial(ch, :).ppc1);
                stsVar(:, :, ii, ch) = std(data_var, 1, 3).*sqrt(size(data_var, 3)-1);
            end
        end
        shaded = 1;
        
    end
    %     txlim = [min(tx) max(tx)];
    txlim = [10 max(tx)];
    tsel = tx>20;
elseif strcmp(thisFile, 'timelockLFP')
    tlabel = allThis(1).label;
    if strcmp(allCfg.name, 'Isis')
        tlabel = cellfun(@ (x) (x(7:end)), tlabel, 'UniformOutput', false);
    elseif strcmp(allCfg.name, 'Ares')
        tlabel = cellfun(@ (x) (x(7:end)), tlabel, 'UniformOutput', false);
    elseif strcmp(allCfg.name, 'Hermes')
        tlabel = cellfun(@ (x) (x(4:end)), tlabel, 'UniformOutput', false);
        tlabel(find(strncmp(tlabel, 'X', 1))) = {'127'};
    end
    tx = allThis(1).time;
    % check if it is concattable
    sizeToCat = unique(cellfun(@(x) (length(x)), {allThis.avg}));
    if length(sizeToCat)>1
        data = cellfun(@(x) (x(:, 1:min(sizeToCat))), {allThis.avg}, 'UniformOutput', false);
        data = cat(3, data{:});
        %         data_var = cellfun(@(x) (x(:, 1:min(sizeToCat))), {allThis.var}, 'UniformOutput', false);
        %         data_var = sqrt(cat(3, data_var{:}));
        tx = tx(1:sizeToCat);
    else
        data = cat(3, allThis.avg);
        %         data_var = sqrt(cat(3, allThis.var));
    end
    %     shaded = 1;
    txlim = [min(tx) max(tx)];
elseif strcmp(thisFile, 'timelockMUAX')
    tlabel = allThis(1).label;
    if strcmp(allCfg.name, 'Isis')
        tlabel = cellfun(@ (x) (x(8:end)), tlabel, 'UniformOutput', false);
    elseif strcmp(allCfg.name, 'Ares')
        tlabel = cellfun(@ (x) (x(8:end)), tlabel, 'UniformOutput', false);
    elseif strcmp(allCfg.name, 'Hermes')
        tlabel = cellfun(@ (x) (x(4:end)), tlabel, 'UniformOutput', false);
        tlabel(find(strncmp(tlabel, 'X', 1))) = {'127'};
    end
    tx = allThis(1).time;
    sizeToCat = unique(cellfun(@(x) (length(x)), {allThis.avg}));
    if length(sizeToCat)>1
        data = cellfun(@(x) (x(:, 1:min(sizeToCat))), {allThis.avg}, 'UniformOutput', false);
        data = cat(3, data{:});
        %         data_var = cellfun(@(x) (x(:, 1:min(sizeToCat))), {allThis.var}, 'UniformOutput', false);
        %         data_var = sqrt(cat(3, data_var{:}));
        tx = tx(1:sizeToCat);
    else
        data = cat(3, allThis.avg);
        %         data_var = sqrt(cat(3, allThis.var));
    end
    %     shaded = 1;
    txlim = [min(tx) max(tx)];
elseif strcmp(thisFile, 'trialPSTH')
    tlabel = allThis(1).label;
    if strcmp(allCfg.name, 'Isis')
        tlabel = cellfun(@ (x) (x(6:end)), tlabel, 'UniformOutput', false);
    elseif strcmp(allCfg.name, 'Ares')
        tlabel = cellfun(@ (x) (x(6:end)), tlabel, 'UniformOutput', false);
    elseif strcmp(allCfg.name, 'Hermes')
        tlabel = cellfun(@ (x) (x(4:end)), tlabel, 'UniformOutput', false);
        tlabel(find(strncmp(tlabel, 'X', 1))) = {'127'};
    end
    tx = allThis(1).time;
    data = cat(3, allThis.avg);
    
    if allCfg.withErrorBars
        data_var = cat(2, allFiles.trialPSTHVar);
        data_var = sqrt(cat(3, data_var.var)/length(data_var));
        shaded = 1;
    end
%     txlim = [min(tx) max(tx)];
txlim = [min(tx) 1.2];
else
    warning('Nothing to plot / not defined')
    return
end

nCond = length(allFiles);
nChan = length(tlabel);

if isfield(allFiles(1), 'sortInd'), sortInd = allFiles(1).sortInd; else sortInd = 1:nCond; end;
assert(length(sortInd)==nCond, 'sortInd doesnt match nCond');
% Channel or Stimuli
if strcmp(allCfg.layout, 'channels')
    nr = size(layout, 1); nc = size(layout, 2);
    for cnd=1:nCond
        cnd = sortInd(cnd);
        h = figure; if allCfg.print; set(h, 'visible', 'off'); end;
        %         him = figure; if allCfg.print; set(him, 'visible', 'off'); end
        for ch=1:nChan
            if caccept(ch)
                thisChan = str2num(tlabel{ch});
                nd = find((layout==thisChan)');
                subplot(nr, nc, nd);
                if strcmp(thisFile, 'stSpec')
                    neighbors = layout(neighbourND(find(layout==thisChan), size(layout), [1 1]));
                    nind = find(ismember(cellfun(@(x) str2num(x), tlabel), neighbors));
                    data = allThis(ch, cnd);
                    data = cat(3, data.ppc1);
                    data_var = stsVar(:, :, cnd, ch);
                    
                    if shaded
                        %                         for ii=nind'
                        %                         shadedErrorBar(tx, data(ii, :), data_var(ii, :)); hold on;
                        %                         end
                        plot(tx, data(nind, :), 'Color', [0.5 0.5 0.5]); hold on;
                        
                    else
                        plot(tx, data(nind, :), 'Color', [0.5 0.5 0.5]); hold on;
                    end
                    %                 ylim([min(min(data([ch nind'], :))) max(max(data([ch nind'], :)))]);
                    if shaded
                        shadedErrorBar(tx, data(ch, :), data_var(ch, :), 'lineprops', 'r', 'transparent', 1); hold on;
                    else
                        plot(tx, data(ch, :), 'r'); hold on;
                    end
                    %                     set(gca, 'XTick', []); set(gca, 'YTick', []);
                    sel = 10 < tx & tx < 120;
                    if allCfg.normalize
                        if shaded
                            ylim([min(min(data(:, sel)-data_var(:, sel))) ...
                                max(max(data(:, sel)+data_var(:, sel)))]);
                        else
                            ylim([min(min(data(:, sel, cnd))) ...
                                max(max(data(:, sel, cnd)))]);
                        end
                    else
                        ylim([min(data(ch, sel, cnd)-data_var(ch, sel, cnd)) ...
                            max(data(ch, sel, cnd)+data_var(ch, sel, cnd))]);
                    end
                    
                else
                    if shaded
                        %                         shadedErrorBar(tx, data(ch, :, cnd), data_var(ch, :, cnd), 'r', 1); hold on;
                        shadedErrorBar(tx, data(ch, :, cnd), data_var(ch, :, cnd), 'lineprops', 'r', 'transparent', 0); hold on;
                        if allCfg.normalize
                            ylim([min(min(data(find(caccept), :, cnd)-data_var(find(caccept), :, cnd))) ...
                                max(max(data(find(caccept), :, cnd)+data_var(find(caccept), :, cnd)))]);
                        else
                            ylim([min(data(ch, sel, cnd)-data_var(ch, :, cnd)) ...
                                max(data(ch, sel, cnd)+data_var(ch, :, cnd))]);
                        end
                    else
                        plot(tx, data(ch, :, cnd), 'r'); hold on;
                        if allCfg.normalize
                            ylim([min(min(data(find(caccept), :, cnd))) ...
                                max(max(data(find(caccept), :, cnd)))]);
                        else
                            ylim([min(data(ch, :, cnd)) max(data(ch, :, cnd))]);
                        end
                    end
                end
                xlim(txlim);
                %                 title(tlabel(ch));
                set(gca,'FontSize',5)
            end
        end
        if allCfg.print
            figname = sprintf('cond%02d_%s.%s', cnd, thisFile, allCfg.ext);
            if strcmp(allCfg.ext, 'eps')
                print(h, fullfile(savefile, figname), sprintf('-d%sc', allCfg.ext), '-r300');
            else
                print(h, fullfile(savefile, figname), sprintf('-d%s', allCfg.ext), '-r300');
            end
        end
        close(h);
    end
else
    if strcmp(allCfg.type, 'grating-ori')
        nr = 6; nc = 12;
    else
        nr = ceil(sqrt(nCond)); nc = ceil(sqrt(nCond));
    end
    for ch=1:nChan
        h = figure(); if allCfg.print; set(h, 'visible', 'off'); end;
        for cnd=1:nCond
            if strcmp(allCfg.type, 'grating-ori')
%                 tok = cellfun(@(x) strsplit(x, '/'), {allFiles.imName}, 'UniformOutput', false);
                tok = cellfun(@(x) strsplit(x, '/'), {allFiles.condName}, 'UniformOutput', false);
                tok = vertcat(tok{:});
                sp = cellfun(@(x) strsplit(x, '_'), tok(:, end), 'UniformOutput', false);
                sp = vertcat(sp{:});
                sori = cellfun(@(x) str2num(x), sp(:, 2), 'UniformOutput', false);
                sfreq = cellfun(@(x) str2num(x(1:3)), sp(:, 3), 'UniformOutput', false);
                orivec = 0:15:165; sfvec = 0.5:0.5:3;
                ndLib = ((length(orivec)*(cellfun(@(x) (find(x == sfvec)), sfreq)-1))...
                    +cellfun(@(x) (find(x == orivec)), sori));
                nd = ndLib(cnd);
            else
                nd = cnd;
            end
            subplot(nr, nc, nd);
            cnd = sortInd(cnd);
            if strcmp(thisFile, 'stSpec')
                thisChan = str2num(tlabel{ch});
                neighbors = layout(neighbourND(find(layout==thisChan), size(layout), [1 1]));
                %                 neighbors = find(neighbors(~isnan(neighbors)));
                nind = find(ismember(cellfun(@(x) str2num(x), tlabel), neighbors));
                data = allThis(ch, cnd);
                datalim = cat(3, allThis(ch, :).ppc1);
                data = cat(3, data.ppc1);
                data_var = stsVar(:, :, cnd, ch);
                datalim_var = stsVar(:, :, :, ch);
                if shaded
                    %                     for ii=nind'
                    %                         shadedErrorBar(tx, data(ii, :), data_var(ii, :)); hold on;
                    %                     end
                    plot(tx, data(nind, :), 'Color', [0.5 0.5 0.5]); hold on;
                    
                else
                    plot(tx, data(nind, :), 'Color', [0.5 0.5 0.5]); hold on;
                end
                
                if shaded
                    shadedErrorBar(tx, data(ch, :), data_var(ch, :), 'lineprops', 'r', 'transparent', 1); hold on;
                else
                    plot(tx, data(ch, :), 'r'); hold on;
                end
                if shaded
                    ylim([min(min(min(datalim([ch nind'], tsel)-datalim_var([ch nind'], tsel)))) ...
                        max(max(max(datalim([ch nind'], tsel)+datalim_var([ch nind'], tsel))))]);
                else
                    ylim([min(min(data([ch nind'], :))) max(max(data([ch nind'], :)))]);
                end
            else
                thisChan = str2num(tlabel{ch});
                if shaded
                    shadedErrorBar(tx, data(ch, :, cnd), data_var(ch, :, cnd), 'lineprops', 'r', 'transparent', 1); hold on;
                    if allCfg.normalize
                        ylim([min(min(data(ch, :, :)-data_var(ch, :, :))) ...
                            max(max(data(ch, :, :)+data_var(ch, :, :)))]);
                    else
                        ylim([min(data(ch, :, cnd)-data_var(ch, :, cnd)) ...
                            max(data(ch, :, cnd)+data_var(ch, :, cnd))]);
                    end
                else
                    plot(tx, data(ch, :, cnd), 'r'); hold on;
                    if allCfg.normalize
                        ylim([min(min(data(ch, :, :))) max(max(data(ch, :, :)))]);
                    else
                        ylim([min(data(ch, :, cnd)) max(data(ch, :, cnd))]);
                    end
                end
                set(gca,'LooseInset',get(gca,'TightInset'))
            end
            xlim(txlim);
            title(sprintf('cond%02d', cnd),'FontSize',5 , 'FontWeight', 'bold');
            set(gca,'FontSize',5)
        end
        if allCfg.print
            figname = sprintf('ch%02d_%s.%s', thisChan, thisFile, allCfg.ext);
            print(h, fullfile(savefile, figname), sprintf('-d%s', allCfg.ext), '-r300');
        end
        close(h);
    end
end
