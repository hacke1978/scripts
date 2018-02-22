function plotLFPpower(allCfg, allFiles, thisFile, layout, varargin)
% plots LFP Power in various forms
addpath('/mnt/v7k/home/uranc/workspace/Source/export_fig/')
savename = allCfg.outputfile;
% flatten data
powerLFP = cat(2, allFiles.(thisFile));

if nargin > 4
    overlayFiles = varargin{1};
    overlayLFP = cat(2, overlayFiles.(thisFile));
    overlayData = cat(3, overlayLFP.powspctrm);
    overlayBaseline = cat(2, overlayFiles.powerBaseline); % DONT FORGET
    overlayBase = cat(3, overlayBaseline.powspctrm);
    overlayBaseline = overlayBase(:, :, 1);
end

nCond = length(powerLFP);
xLab = powerLFP.freq;
yLab = powerLFP.label;
nChan = length(yLab);
data = cat(3, powerLFP.powspctrm);

if isfield(allFiles(1), 'sortInd'), sortInd = allFiles(1).sortInd; else sortInd = 1:nCond; end;
assert(length(sortInd)==nCond, 'sortInd doesnt match nCond');

% load if possible
if isfield(allFiles, 'lfpPower2')
    powerLFP2 = cat(2, allFiles.lfpPower2);
    data2 = cat(3, powerLFP2.powspctrm);
    data = (data + data2)/2;
    disp('averaging two windows!!!');
end

% load if possible
if isfield(allFiles, 'lfpPower_Second')
    powerSec = cat(2, allFiles.lfpPower_Second);
    data_Second = cat(3, powerSec.powspctrm);
end
caccept = find(allFiles(1).caccept);
% get baseline
powerBaseline = cat(2, allFiles.powerBaseline); % DONT FORGET
base = cat(3, powerBaseline.powspctrm);
% powerBaseline = cat(2, allFiles(end).(thisFile));
% baseline = cat(3, powerBaseline.powspctrm);

if allCfg.withErrorBars
    lfpVar = zeros([size(allFiles(1).lfpPowerVar(1).powspctrm) nCond]);
    for ii=1:nCond
        if isfield(allCfg, 'nConds')
            baseline = base(:, :, sum(ii>cumsum(allCfg.nConds))+1);
        else
            baseline = base(:, :, 1);
        end
        tVar = cat(3, allFiles(ii).lfpPowerVar.powspctrm);
        tVar = log10(tVar./repmat(baseline, 1, 1, size(tVar, 3)));
        lfpVar(:,:,ii) = std(tVar, 1, 3).*sqrt(size(tVar, 3)-1);
    end
end

allNormalized = nan(size(data(:, :, 1)));
for cnd =1:nCond
    if isfield(allCfg, 'nConds')
        baseline = base(:, :, sum(cnd>cumsum(allCfg.nConds))+1);
    else
        baseline = base(:, :, 1);
    end
    allNormalized = [allNormalized mean(data(:, :, cnd)./baseline, 2)];
    
end
[~, b] = min(allNormalized, [], 2); probGrayInd = mode(b); %% Dangerous but logical

% Dimensions for the subplot
if strcmp(allCfg.layout, 'channels')
    nr = size(layout, 1); nc = size(layout, 2);
    for cnd=1:nCond
        if isfield(allCfg, 'nConds')
            baseline = base(:, :, sum(cnd>cumsum(allCfg.nConds))+1);
        else
            baseline = base(:, :, 1);
        end
        cnd = sortInd(cnd);
        h1 = figure('Visible', 'off'); %if allCfg.print; set(h1, ); end;
        h2 = figure('Visible', 'off'); %if allCfg.print; set(h2, 'Visible', 'off'); end;
        h3 = figure('Visible', 'off'); %if allCfg.print; set(h3, 'Visible', 'off'); end;
        %         if allCfg.onIm; him = figure; if allCfg.print; set(him, 'visible', 'off'); end; end;
        for ch=1:nChan
            if strcmp(allCfg.name, 'Isis')
                chn = str2num(yLab{ch}(7:end));
            elseif strcmp(allCfg.name, 'Ares')
                chn = str2num(yLab{ch}(7:end));
            elseif strcmp(allCfg.name, 'Hermes')
                chn = str2num(yLab{ch}(4:end));
                if isempty(chn)
                    chn = 127;
                end
            end
            if chn == 84 || chn == 127
                continue;
            end
            nd = find((layout==chn)');
            %raw plot
            %             figure(h1); if allCfg.print; set(h1, 'Visible', 'off'); end;
            set(0,'CurrentFigure',h1);
            subplot(nr, nc, nd);
            plot(xLab, baseline(ch, :), 'Color', [0.5 0.5 0.5]);
            plot(xLab, data(ch, :, cnd), 'r'); hold on;
            title({chn}, 'FontWeight', 'bold', 'FontSize', 5);
            ylim([min([data(ch, :, cnd) baseline(ch, :)]) max([data(ch, :, cnd) baseline(ch, :)])]);
            xlim([10 140]); %ylim([-0.5 1.2]);
            set(gca,'FontSize',5);
            
            % baseline normalized
            %             figure(h2); if allCfg.print; set(h2, 'Visible', 'off'); end;
            set(0,'CurrentFigure',h2);
            subplot(nr, nc, nd);
            if allCfg.withErrorBars
                shadedErrorBar(xLab, log10(data(ch, :, cnd)./(1e-20+baseline(ch, :))), ...
                    lfpVar(ch, :, cnd), 'lineprops', 'r', 'transparent', 1); hold on;
            else
                plot(xLab, log10(data(ch, :, cnd)./(1e-20+baseline(ch, :))), 'r'); hold on;
            end
            %             if allCfg.onIm; plotOnIm(him, xLab, log10(data(ch, :, cnd)./baseline(ch, :)); end;
            title({chn}, 'FontWeight', 'bold', 'FontSize', 5);
            xlim([10 140]); %ylim([-0.2 1]);
            sel = 10 < xLab & xLab < 140;
            if allCfg.normalize
                if allCfg.withErrorBars
                    ylim([min(min(log10((data(caccept, sel, cnd))./(1e-20+baseline(caccept, sel)))-lfpVar(caccept, sel, cnd))) ...
                        max(max(log10((data(caccept, sel, cnd))./(1e-20+baseline(caccept, sel)))+lfpVar(caccept, sel, cnd)))]);
                else
                    ylim([min(min(log10(data(caccept, sel, cnd)./(1e-20+baseline(caccept, sel))))) ...
                        max(max(log10(data(caccept, sel, cnd)./(1e-20+baseline(caccept, sel)))))]);
                end
            else
                ylim([min(log10(data(ch, sel, cnd)./(1e-20+baseline(ch, sel)))) ...
                    max(log10(data(ch, sel, cnd)./(1e-20+baseline(ch, sel))))]);
            end
            set(gca,'FontSize',5)
            %             if ~strcmp(allCfg.type, 'grating-ori')
            %                 % gray normalized
            %                 figure(h3); if allCfg.print; set(h3, 'Visible', 'off'); end;
            %                 subplot(nr, nc, nd);
            %                 if allCfg.withErrorBars
            %                     shadedErrorBar(xLab, log10(data(ch, :, cnd)./(1e-20+baseline(ch, :))), lfpVar(ch, :, cnd), 'r', 1); hold on;
            %                 else
            %                     plot(xLab, log10(data(ch, :, cnd)./(1e-20+data(ch, :, probGrayInd))), 'r'); hold on;
            %                 end
            %                 %             ylim([min(log10(data(ch, :, cnd)./(1e-20+data(ch, :, probGrayInd)))) max(log10(data(ch, :, cnd)./(1e-20+data(ch, :, probGrayInd))))]);
            %                 sel = 10 < xLab & xLab < 120;
            %                 if allCfg.normalize
            %                     if allCfg.withErrorBars
            %                         ylim([min(min((log10(data(caccept, sel, cnd)./(data(caccept, sel, probGrayInd)))-lfpVar(caccept,sel,cnd)))) ...
            %                             max(max(log10((data(caccept, sel, cnd)./(1e-20+data(caccept, sel, probGrayInd))))+lfpVar(caccept,sel,cnd)))]);
            %                     else
            %                         ylim([min(min(log10(data(caccept, sel, cnd)./(1e-20+data(caccept, sel, probGrayInd))))) ...
            %                             max(max(log10(data(caccept, sel, cnd)./(1e-20+data(caccept, sel, probGrayInd)))))]);
            %                     end
            %                 else
            %                     ylim([min(log10(data(ch, :, cnd)./(1e-20+baseline(ch, :)))) ...
            %                         max(log10(data(ch, :, cnd)./(1e-20+baseline(ch, :))))]);
            %                 end
            % %                 legend({num2str(chn)})
            % %                 title({chn}, 'FontWeight', 'bold', 'FontSize', 5);
            %                 %ylim([-0.2 1]);
            %                 xlim([10 120]);
            %                 set(gca,'FontSize',5);
            %             end
        end
        
        % print all
        if allCfg.print
            if isfield(allFiles, 'lfpPower2')
                fname = sprintf('cond%02d_%s_raw2.%s', cnd, thisFile, allCfg.ext);
            else
                fname = sprintf('cond%02d_%s_raw.%s', cnd, thisFile, allCfg.ext);
            end
            print(h1, fullfile(savename, fname), sprintf('-d%s', allCfg.ext), '-r300');
            if isfield(allFiles, 'lfpPower2')
                fname = sprintf('cond%02d_%s_base_normalized2.%s', cnd, thisFile, allCfg.ext);
            else
                fname = sprintf('cond%02d_%s_base_normalized.%s', cnd, thisFile, allCfg.ext);
            end
            print(h2, fullfile(savename, fname), sprintf('-d%s', allCfg.ext), '-r300');
            %             if ~strcmp(allCfg.type, 'grating-ori')
            %                 if isfield(allFiles, 'lfpPower2')
            %                     fname = sprintf('cond%02d_%s_gray_normalized2.%s', cnd, thisFile, allCfg.ext);
            %                 else
            %                     fname = sprintf('cond%02d_%s_gray_normalized.%s', cnd, thisFile, allCfg.ext);
            %                 end
            %                 print(h3, fullfile(savename, fname), sprintf('-d%s', allCfg.ext), '-r300');
            %             end
        end
        close(h1); close(h2); close(h3);
    end
elseif strcmp(allCfg.layout, 'stimuli')
    if strcmp(allCfg.type, 'grating-ori')
        nr = 6+1; nc = 12+1;
    else
        nr = ceil(sqrt(nCond)); nc = ceil(sqrt(nCond));
    end
    for ch=1:nChan
        if strcmp(allCfg.name, 'Isis')
            chn = str2num(yLab{ch}(7:end));
        elseif strcmp(allCfg.name, 'Ares')
            chn = str2num(yLab{ch}(5:end));
        elseif strcmp(allCfg.name, 'Hermes')
            chn = str2num(yLab{ch}(4:end));
            if isempty(chn)
                chn = 127;
            end
        end
        h1 = figure('Visible', 'off');% if allCfg.print; set(h1, 'Visible', 'off'); end;
        h2 = figure('Visible', 'off');% if allCfg.print; set(h2, 'Visible', 'off'); end;
        h3 = figure('Visible', 'off');% if allCfg.print; set(h3, 'Visible', 'off'); end;
        for cnd=1:nCond
            if isfield(allCfg, 'nConds')
                baseline = base(:, :, sum(cnd>cumsum(allCfg.nConds))+1);
            else
                baseline = base(:, :, 1);
            end
            if strcmp(allCfg.type, 'grating-ori')
                ndLayout = reshape([1:nr*nc], nc, nr)';
                %                 tok = cellfun(@(x) strsplit(x, '/'), {allFiles.imName}, 'UniformOutput', false);
                tok = cellfun(@(x) strsplit(x, '/'), {allFiles.condName}, 'UniformOutput', false);
                tok = vertcat(tok{:});
                sp = cellfun(@(x) strsplit(x, '_'), tok(:, end), 'UniformOutput', false);
                sp = vertcat(sp{:});
                sori = cellfun(@(x) str2num(x), sp(:, 2), 'UniformOutput', false);
                sfreq = cellfun(@(x) str2num(x(1:3)), sp(:, 3), 'UniformOutput', false);
                orivec = 0:15:165; sfvec = 0.5:0.5:3;
                ndLib = (((length(orivec)+1)*(cellfun(@(x) (find(x == sfvec)), sfreq)-1))...
                    +cellfun(@(x) (find(x == orivec)), sori));
                nd = ndLib(cnd);
            else
                nd = cnd;
            end
            cnd = sortInd(cnd);
            % raw
            %             figure(h1); if allCfg.print; set(h1, 'Visible', 'off'); end;
            set(0,'CurrentFigure',h1);
            subplot(nr, nc, nd)
            %             plot(xLab, baseline(ch, :), 'Color', [0.5 0.5 0.5]); hold on;
            if allCfg.isOverlay
                plot(xLab, overlayBaseline(ch, :), 'Color', [0 0 0.5]); hold on;
            end
            loglog(xLab, (baseline(ch, :)), 'Color', [0.5 0.5 0.5], 'linewidth', 1); hold on;
            if strcmp(allCfg.type, 'NatImSEQ')
                if nd<nCond/2
                    plot(xLab, data(ch, :, cnd), 'k'); hold on;
                    plot(xLab, data_Second(ch, :, cnd), 'r'); hold on;
                else
                    plot(xLab, data(ch, :, cnd), 'r'); hold on;
                    plot(xLab, data_Second(ch, :, cnd), 'k'); hold on;
                end
            else
                %                 plot(xLab, data(ch, :, cnd), 'r'); hold on;
                if allCfg.isOverlay
                    plot(xLab, overlayData(ch, :, cnd), 'b'); hold on;
                end
                loglog(xLab, (data(ch, :, cnd)), 'r', 'linewidth', 1); hold on;
            end
            xlim([10 140]);
            if allCfg.gammaPeak
                if strcmp(allCfg.gammaPeak, 'all');
                    [base_bias, base_exp, event_bias, event_exp_2, gauss_amp, gauss_freq, gauss_std, fit_f2] = ...
                        fit_gammadata(round(xLab)', round(xLab), base(ch, :), data(ch, :, cnd));
                else
%                     [base_exp, base_bias,event_bias, event_exp,exp_coeff, gauss_amp, gauss_freq, fit_f2] = ...
%                         fit_gammadata(round(xLab)', allCfg.gammaPeak, base(ch, :), data(ch, :, cnd));
                [base_bias, base_exp, event_bias, event_exp_2, gauss_amp, gauss_freq, gauss_std, fit_f2] = ...
                        fit_gammadata(round(xLab)', allCfg.gammaPeak, base(ch, :), data(ch, :, cnd));
                end                
%                 figure
%                 loglog((xLab), data(ch, :, cnd), 'color', [1 0 0 ]/2), hold on
             %   keyboard
                loglog((xLab), 10.^(base_bias-log10(round(xLab))*base_exp), 'color', [.5 .5 .5]/2)
                loglog((xLab), 10.^(fit_f2), 'color', [1 0 0 ]/2), hold on
%                 loglog((xLab), 10.^(fit_linear), 'color', [1 0 0 ]), hold on
%                 loglog((xLab), 10.^(fit_gauss_1), 'color', [0 1 0 ]), hold on
%                 loglog((xLab), 10.^(fit_gauss_2), 'color', [0 0 1 ]), hold on
%                 loglog((xLab), 10.^(event_bias-log10(round(xLab))*event_exp+exp(-log10(round(xLab))*exp_coeff)), 'color', [1 0 0 ]), hold on 
                text((10^gauss_freq)/2, min(10.^(fit_f2))/10,...
                    sprintf('Peak: %2.2f\nFreq: %2.2f\nstd: %2.2f', gauss_amp, 10^gauss_freq, 10^gauss_std),...
                    'fontsize', 5); hold on;
            end
            
            %             ylim([min([data(ch, :, cnd) baseline(ch, :)]) max([data(ch, :, cnd) baseline(ch, :)])]);
            if strcmp(allCfg.type, 'grating-ori')
                if mod(ndLib(cnd), nc)==1
                    title(sprintf('%.2f', sfreq{cnd}), 'FontWeight', 'bold', 'FontSize', 3);
                elseif ceil(ndLib(cnd)/nc)==1
                    title(sprintf('%d', sori{cnd}), 'FontWeight', 'bold', 'FontSize', 3);
                elseif allCfg.normalize
                    set(gca, 'XTick', []); set(gca, 'YTick', []);
                end
            end
            set(gca,'FontSize', 3)
            
            % normalized to baseline
            %             figure(h2); if allCfg.print; set(h2, 'Visible', 'off'); end;
            set(0,'CurrentFigure',h2);
            subplot(nr, nc, nd)
            plot(xLab, mean(squeeze(log10(data(ch, :, :)./repmat(1e-20+baseline(ch, :), 1, 1, nCond))), 2), '--', 'Color', [0.5 0.5 0.5]); hold on;
            if allCfg.withErrorBars
                shadedErrorBar(xLab, log10(data(ch, :, cnd)./(1e-20+baseline(ch, :))), ...
                    lfpVar(ch, :, cnd), 'lineprops', 'r', 'transparent', 1); hold on;
            else
                if strcmp(allCfg.type, 'NatImSEQ')
                    if nd<nCond/2
                        plot(xLab, log10(data(ch, :, cnd)./(1e-20+baseline(ch, :))), 'k'); hold on;
                        plot(xLab, log10(data_Second(ch, :, cnd)./(1e-20+baseline(ch, :))), 'r'); hold on;
                    else
                        plot(xLab, log10(data(ch, :, cnd)./(1e-20+baseline(ch, :))), 'r'); hold on;
                        plot(xLab, log10(data_Second(ch, :, cnd)./(1e-20+baseline(ch, :))), 'k'); hold on;
                    end
                else
                    if allCfg.isOverlay
                        plot(xLab, log10(overlayData(ch, :, cnd)./(1e-20+overlayBase(ch, :))), 'b'); hold on;
                    end
                    plot(xLab, log10(data(ch, :, cnd)./(1e-20+baseline(ch, :))), 'r'); hold on;
                end
            end
            %             title(sprintf('cond%02d', cnd), 'FontWeight', 'bold', 'FontSize', 5);
            if strcmp(allCfg.type, 'grating-ori')
                if mod(ndLib(cnd), nc)==1
                    title(sprintf('%.2f', sfreq{cnd}), 'FontWeight', 'bold', 'FontSize', 3);
                elseif ceil(ndLib(cnd)/nc)==1
                    title(sprintf('%d', sori{cnd}), 'FontWeight', 'bold', 'FontSize', 3);
                elseif allCfg.normalize
                    set(gca, 'XTick', []); set(gca, 'YTick', []);
                end
            end
            xlim([10 140]);
            set(gca,'FontSize', 3)
            sel = 10 < xLab & xLab < 140;
            if allCfg.normalize
                ylim([min(min(log10(data(ch, sel, :)./repmat(1e-20+baseline(ch, sel), 1, 1, size(data, 3))))) ...
                    max(max(log10(data(ch, sel, :)./repmat(1e-20+baseline(ch, sel), 1, 1, size(data, 3)))))]);
            else
                ylim([min(log10(data(ch, sel, cnd)./(1e-20+baseline(ch, sel)))) ...
                    max(log10(data(ch, sel, cnd)./(1e-20+baseline(ch, sel))))]);
            end
            %             ylim([min(log10(data(ch, :, cnd)./(1e-20+baseline(ch, :)))) max(log10(data(ch, :, cnd)./(1e-20+baseline(ch, :))))])
            
            % normalized to gray
            %             if ~strcmp(allCfg.type, 'grating-ori')
            %                 figure(h3); if allCfg.print; set(h3, 'Visible', 'off'); end;
            %                 subplot(nr, nc, nd)
            %                 if allCfg.withErrorBars
            %                     shadedErrorBar(xLab, log10(data(ch, :, cnd)./(1e-20+data(ch, :, probGrayInd))), ...
            %                         lfpVar(ch, :, cnd), 'lineprops', 'r', 'transparent', 1); hold on;
            %                 else
            %                     plot(xLab, log10(data(ch, :, cnd)./(1e-20+data(ch, :, probGrayInd))), 'r'); hold on;
            %                 end
            %                 title(sprintf('cond%02d', cnd), 'FontWeight', 'bold', 'FontSize', 5);
            %                 xlim([10 120]); %ylim([-0.2 1]);
            %                 if allCfg.normalize
            %                     if allCfg.withErrorBars
            %                         ylim([min(min(log10((data(ch, sel, :))./repmat(1e-20+baseline(ch, sel), 1, 1, size(data, 3)))-lfpVar(ch,sel,:))) ...
            %                             max(max(log10((data(ch, sel, :))./repmat(1e-20+baseline(ch, sel), 1, 1, size(data, 3)))+lfpVar(ch,sel,:)))]);
            %                     else
            %                         ylim([min(min(log10(data(ch, sel, :)./repmat(1e-20+baseline(ch, sel), 1, 1, size(data, 3))))) ...
            %                             max(max(log10(data(ch, sel, :)./repmat(1e-20+baseline(ch, sel), 1, 1, size(data, 3)))))]);
            %                     end
            %                 else
            %                     ylim([min(log10(data(ch, sel, cnd)./(1e-20+baseline(ch, sel)))) ...
            %                         max(log10(data(ch, sel, cnd)./(1e-20+baseline(ch, sel))))]);
            %                 end
            %                 %             ylim([min(log10(data(ch, :, cnd)./(1e-20+data(ch, :, probGrayInd)))) max(log10(data(ch, :, cnd)./(1e-20+data(ch, :, probGrayInd))))]);
            %                 set(gca, 'XTick', []); set(gca, 'YTick', []);
            %                 set(gca,'FontSize', 5)
            %             end
        end
        if strcmp(allCfg.type, 'grating-ori')
            ndLayout = reshape([1:nr*nc], nc, nr)';
            for rr=1:nr-1
                subplot(nr, nc, rr*nc)
                plot(xLab, mean(squeeze(log10(data(ch, :, :)./repmat(1e-20+baseline(ch, :), 1, 1, nCond))), 2), '--', 'Color', [0.5 0.5 0.5]); hold on;
                cind = cellfun(@(x) find(x==ndLib), num2cell(ndLayout(rr, 1:nc-1)), 'UniformOutput', false);
                plot(xLab, log10(mean(data(ch, :, [cind{:}]), 3)./(1e-20+baseline(ch, :))), 'r'); hold on;
                if allCfg.isOverlay
                    plot(xLab, log10(mean(overlayData(ch, :, [cind{:}]), 3)./(1e-20+overlayBase(ch, :))), 'b'); hold on;
                end
                xlim([10 140]); %ylim([-0.2 1]);
                set(gca, 'XTick', []); set(gca, 'YTick', []);
                set(gca,'FontSize', 3)
                
                sel = 10 < xLab & xLab < 140;
                if allCfg.normalize
                    ylim([min(min(log10(data(ch, sel, :)./repmat(1e-20+baseline(ch, sel), 1, 1, size(data, 3))))) ...
                        max(max(log10(data(ch, sel, :)./repmat(1e-20+baseline(ch, sel), 1, 1, size(data, 3)))))]);
                else
                    ylim([min(log10(data(ch, sel, cnd)./(1e-20+baseline(ch, sel)))) ...
                        max(log10(data(ch, sel, cnd)./(1e-20+baseline(ch, sel))))]);
                end
            end
            for cc=1:nc-1
                subplot(nr, nc, nc*(nr-1)+cc)
                plot(xLab, mean(squeeze(log10(data(ch, :, :)./repmat(1e-20+baseline(ch, :), 1, 1, nCond))), 2), '--', 'Color', [0.5 0.5 0.5]); hold on;
                cind = cellfun(@(x) find(x==ndLib), num2cell(ndLayout(1:nr-1, cc)), 'UniformOutput', false);
                plot(xLab, log10(mean(data(ch, :, [cind{:}]), 3)./(1e-20+baseline(ch, :))), 'r'); hold on;
                if allCfg.isOverlay
                    plot(xLab, log10(mean(overlayData(ch, :, [cind{:}]), 3)./(1e-20+overlayBase(ch, :))), 'b'); hold on;
                end
                xlim([10 140]); %ylim([-0.2 1]);
                set(gca, 'XTick', []); set(gca, 'YTick', []);
                set(gca,'FontSize', 3)
                
                sel = 10 < xLab & xLab < 140;
                if allCfg.normalize
                    ylim([min(min(log10(data(ch, sel, :)./repmat(1e-20+baseline(ch, sel), 1, 1, size(data, 3))))) ...
                        max(max(log10(data(ch, sel, :)./repmat(1e-20+baseline(ch, sel), 1, 1, size(data, 3)))))]);
                else
                    ylim([min(log10(data(ch, sel, cnd)./(1e-20+baseline(ch, sel)))) ...
                        max(log10(data(ch, sel, cnd)./(1e-20+baseline(ch, sel))))]);
                end
            end
        end
        % print all
        if allCfg.print
            if isfield(allFiles, 'lfpPower2')
                fname = sprintf('ch%02d_%s_raw2.%s', chn, thisFile, allCfg.ext);
            else
                fname = sprintf('ch%02d_%s_raw.%s', chn, thisFile, allCfg.ext);
            end
            
            print(h1, fullfile(savename, fname), sprintf('-d%s', allCfg.ext), '-r300');
            if isfield(allFiles, 'lfpPower2')
                fname = sprintf('ch%02d_%s_base_normalized2.%s', chn, thisFile, allCfg.ext);
            else
                fname = sprintf('ch%02d_%s_base_normalized.%s', chn, thisFile, allCfg.ext);
            end
            
            print(h2, fullfile(savename, fname), sprintf('-d%s', allCfg.ext), '-r300');
            %             if ~strcmp(allCfg.type, 'grating-ori')
            %                 if isfield(allFiles, 'lfpPower2')
            %                     fname = sprintf('ch%02d_%s_gray_normalized2.%s', chn, thisFile, allCfg.ext);
            %                 else
            %                     fname = sprintf('ch%02d_%s_gray_normalized.%s', chn, thisFile, allCfg.ext);
            %                 end
            %                 print(h3, fullfile(savename, fname), sprintf('-d%s', allCfg.ext), '-r300');
            %             end
        end
        close(h1); close(h2); close(h3);
    end
end
