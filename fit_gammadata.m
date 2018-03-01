function [fit_params, fit_line] = fit_gammadata(f, f_use4fit, data_base, data_fit, fitType)
% input:
% f: frequencies
% f_use4fit=[35:57 65:115 126:175 186:200]; can exclude freq
% data_base: used to fit exp: (1/f^exp) - enter power (not log)
% data_fit: used to fit weights and gaussian - enter power (not log)
% output (fit_params fit_line)

%     Copyright (C) 2014  D Hermes
%
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
%
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
%
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.


if size(f, 2) < size(f, 1); f = f'; end % function expects row vectors
if size(f_use4fit, 2) < size(f_use4fit, 1); f = f'; end

f_sel=ismember(f,f_use4fit);
x_base=data_base(f_sel);
x_in=data_fit(f_sel);
f_in=f(f_sel);

if ~isequal(length(f_in), length(x_base)); x_base = x_base'; end
if ~isequal(length(f_in), length(x_in)); x_in = x_in'; end

% fit exponent to base
p=polyfit(log10(f_in),log10(x_base),1);
base_exp  = -p(1);
base_bias =  p(2);

% keyboard
my_options=optimset('Display','off','Algorithm','trust-region-reflective');

switch fitType
    case 'exp'
        % fit params
        flagStep = true;
        flagFine = true;
        
        % exp params
        f_in_1=f_in; f_in_1(f_in>25 & f_in<120) = [];
        x_in_1 = x_in; x_in_1(f_in>25 & f_in<120) = [];
        
        % first fit
        slopeinit = linspace(1,10,3);
        tauinit = logspace(1e-2,100,30);
        rms = inf;
        x = [];
        if flagStep
            for slope_init = 1:length(slopeinit)
                for tau_init = 1:length(tauinit)
                    [x_tmp, ~, residual]=lsqnonlin(@(x) fit_loglog_exp(x, log10(x_in_1), log10(f_in_1), flagStep),...
                        [ rand(1,1)  -slopeinit(slope_init)     tauinit(tau_init) ],...
                        [ -Inf       -Inf                       0                 ],...
                        [ +Inf        0                         Inf               ],...
                        my_options);
                    if sum(residual.^2)<rms
                        rms = sum(residual.^2);
                        x = x_tmp;
                    end
                end
            end
            exp_param = x;
            
            % fit the gaussians on top of the exponential
            x = [];
            rms = Inf;
            finit = linspace(30,80,5);
            sdinit = linspace(log10(1./sqrt(2/3)), log10(1./sqrt(19/20)),3);
            for iInit = 1:length(finit)
                for jSd_init = 1:length(sdinit)
                    [x_tmp, ~, residual]=lsqnonlin(@(x) fit_loglog_exp(x, log10(x_in), log10(f_in), flagStep, exp_param),...
                        [    0    log10(finit(iInit)) sdinit(jSd_init)       0      log10(1)      ],...
                        [    0    log10(30)           log10(1./sqrt(19/20))  -pi    log10(0.95)   ],...
                        [    Inf  log10(80)           log10(1./sqrt(2/3))    pi     log10(1.1)    ],...
                        my_options);
                    if sum(residual.^2)<rms
                        rms = sum(residual.^2);
                        x = x_tmp;
                    end
                end
            end
            % second gaussian is constrained by the first
            x(6) = x(3); % std
            x(4) = (sin(x(4))+1)*0.5*x(1);
            %             x(4) = x(1)-x(4); % amp
            x(5) = x(2)+log10(2)+x(5); % peak
            
            % % %
            x = [exp_param x;]; % pack all params
            
            %             keyboard
            if flagFine
                % get lines
                fit_linear = (x(1)-x(2)*exp(-log10(f)*x(3)));
                fit_gauss_1 = x(4)*x(6)*sqrt(2*pi)*normpdf(log10(f),x(5),x(6));
                fit_gauss_2 =  x(7)*x(9)*sqrt(2*pi)*normpdf(log10(f),x(8),x(9));
                fit_poly = fit_linear + fit_gauss_1 + fit_gauss_2;
                
                % exp params
                f_in_1=f_in; % second(fine) step includs gamma range[30-100]
                x_in_1 = 10.^(log10(x_in)-(fit_gauss_1(f_sel)+fit_gauss_2(f_sel)));
                
                % first fit
                slopeinit = linspace(1,10,3);
                tauinit = logspace(1e-2,100,30);
                rms = inf;
                x = [];
                
                for slope_init = 1:length(slopeinit)
                    for tau_init = 1:length(tauinit)
                        [x_tmp, ~, residual]=lsqnonlin(@(x) ...
                            fit_loglog_exp(x, log10(x_in_1), log10(f_in_1), flagStep),...
                            [ rand(1,1)  -slopeinit(slope_init)     tauinit(tau_init) ],...
                            [ -Inf       -Inf                       0                 ],...
                            [ +Inf        0                         Inf               ],...
                            my_options);
                        if sum(residual.^2)<rms
                            rms = sum(residual.^2);
                            x = x_tmp;
                        end
                    end
                end
                exp_param = x;
                
                % fit the gaussians on top of the exponential
                x = [];
                rms = Inf;
                finit = linspace(30,80,5);
                sdinit = linspace(log10(1./sqrt(2/3)), log10(1./sqrt(19/20)),3);
                for iInit = 1:length(finit)
                    for jSd_init = 1:length(sdinit)
                        [x_tmp, ~, residual]=lsqnonlin(@(x) ...
                            fit_loglog_exp(x, log10(x_in), log10(f_in), flagStep, exp_param),...
                            [    0    log10(finit(iInit)) sdinit(jSd_init)       0      log10(1)      ],...
                            [    0    log10(30)           log10(1./sqrt(19/20))  -pi    log10(0.95)   ],...
                            [    Inf  log10(80)           log10(1./sqrt(2/3))    pi     log10(1.1)    ],...
                            my_options);
                        if sum(residual.^2)<rms
                            rms = sum(residual.^2);
                            x = x_tmp;
                        end
                    end
                end
                % second gaussian is constrained by the first
                x(6) = x(3); % std
                x(4) = (sin(x(4))+1)*0.5*x(1);
                %             x(4) = x(1)-x(4); % amp
                x(5) = x(2)+log10(2)+x(5); % peak
                
                % % %
                x = [exp_param x;]; % pack all params
            end
        else
            x = [];
            rms = inf;
            finit = linspace(30,80,5);
            sdinit = linspace(log10(1./sqrt(2/3)), log10(1./sqrt(19/20)),3);
            tauinit = logspace(1e-2,100,3);
            for iInit = 1:length(finit)
                for jSd_init = 1:length(sdinit)
                    for kSlope_init = 1:length(tauinit)
                        [x_tmp, ~, residual]=lsqnonlin(@(x) fit_loglog_exp(x, log10(x_in), log10(f_in), false),...
                            [ 0          0    tauinit(kSlope_init)    0   log10(finit(iInit)) sdinit(jSd_init)       0    log10(finit(iInit))+log(2) sdinit(jSd_init)      ],...
                            [ -Inf      -Inf  0                       0   log10(30)           log10(1./sqrt(19/20))  0    log10(60)                  log10(1./sqrt(19/20)) ],...
                            [ +Inf      +Inf  Inf                     Inf log10(80)           log10(1./sqrt(2/3))    Inf  log10(f(end))              log10(1./sqrt(2/3))   ],...
                            my_options);
                        if sum(residual.^2)<rms
                            rms = sum(residual.^2);
                            x = x_tmp;
                        end
                    end
                end
            end
        end
        
        % get params
        fit_params = [];
        fit_params.base_bias = base_bias;
        fit_params.base_exp = base_exp;
        fit_params.stim_bias=x(1);
        fit_params.stim_exp=x(2);
        fit_params.stim_exp_2=x(3);
        fit_params.gauss_amp=x(4);
        fit_params.gauss_freq=x(5);
        fit_params.gauss_std=x(6);
        fit_params.gauss_amp_2=x(7);
        fit_params.gauss_freq_2=x(8);
        fit_params.gauss_std_2=x(9);
        
        % get lines
        fit_linear = (x(1)-x(2)*exp(-log10(f)*x(3)));
        fit_gauss_1 = x(4)*x(6)*sqrt(2*pi)*normpdf(log10(f),x(5),x(6));
        fit_gauss_2 =  x(7)*x(9)*sqrt(2*pi)*normpdf(log10(f),x(8),x(9));
        fit_line = fit_linear + fit_gauss_1 + fit_gauss_2;
        
    case 'linear'
        % fit params
        f_in_1=f_in; f_in_1(f_in_1<40) = nan;
        
        x = [];
        rms = inf;
        finit = linspace(30,80,5);
        sdinit = linspace(log10(1./sqrt(2/3)), log10(1./sqrt(19/20)),3);
        slopeinit = [0 0.25*pi 0.5*pi];
        for iInit = 1:length(finit)
            for jSd_init = 1:length(sdinit)
                for kSlope_init = 1:length(slopeinit)
                    [x_tmp, ~, residual]=lsqnonlin(@(x) fit_func3_loglog(x, log10(x_in), log10(f_in), log10(f_in_1)),...
                        [ base_bias  slopeinit(kSlope_init) 0   log10(finit(iInit)) sdinit(jSd_init)       base_exp 0   log10(1)     ],...
                        [ -Inf       -pi                    0   log10(30)           log10(1./sqrt(19/20))  0        0   log10(0.95)  ],...
                        [ +Inf        pi                    Inf log10(80)           log10(1./sqrt(2/3))    Inf      Inf log10(1.1)   ],...
                        my_options);
                    
                    if sum(residual.^2)<rms
                        rms = sum(residual.^2);
                        x = x_tmp;
                    end
                    
                end
            end
        end
        x(2) = -(sin(x(2))+1)*0.5*x(6);
        x(7) = x(3)-x(7);
        x(8) = x(4)+log10(2)+x(8);
        x(9) = x(5);
        fit_params = [];
        fit_params.base_bias=base_bias;
        fit_params.base_exp=base_exp;
        fit_params.stim_bias=x(1);
        fit_params.stim_exp=x(2);
        fit_params.stim_exp_2=x(6);
        fit_params.gauss_amp=x(3);
        fit_params.gauss_freq=x(4);
        fit_params.gauss_std=x(5);
        fit_params.gauss_amp_2=x(7);
        fit_params.gauss_freq_2=x(8);
        fit_params.gauss_std_2=x(9);
        
        % get line on the whole range
        f_1 = f; f_1(f_1<40) = nan;
        fit_linear = x(1)-nansum([x(6)*log10(f);x(2)*log10(f_1)-x(2)*log10(40)], 1);
        fit_gauss_1 = x(3)*x(5)*sqrt(2*pi)*normpdf(log10(f),x(4), x(5));
        fit_gauss_2 = x(7)*x(9)*sqrt(2*pi)*normpdf(log10(f),x(8), x(9));
        fit_poly = fit_linear + fit_gauss_1 + fit_gauss_2;
    case 'poly'
        
        % fit exponent to base
        peak_bw = 20; % hz
        f_peak = f_in(30 < f_in & f_in < 120); % where to look for peaks
        
        % get best polyfit
        rmsFit = Inf;
        bestP = [];
        for ii=1:30
            [p, S] = polyfit(log10(f_in),log10(x_in), ii);
            if S.normr < rmsFit
                rmsFit = S.normr;
                bestP = p;
                bestOrder = ii;
            end
        end
        fit_poly = polyval(p, log10(f_in));
        
        % find peaks
        peak_min = []; peak_max = [];
        peak_min = peakseek(fit_poly, find(cumsum(diff(f_in))==peak_bw/2), min(fit_poly),'right', 'min');
        peak_max = peakseek(fit_poly, find(cumsum(diff(f_in))==peak_bw/2), min(fit_poly),'right', 'max');
        tdel_max = peak_max < peak_min(1) | f_in(peak_max) < min(f_peak)...
                                          | f_in(peak_max) > max(f_peak);
        tdel_min = peak_min > peak_max(end) | peak_min > max(f_peak);
        peak_max(tdel_max) = [];
        peak_min(tdel_min) = [];
%         if length(peak_min)~=length(peak_max)
            new_max = [];
            for ii=1:length(peak_min)
                if isempty(peak_max(find((peak_max- peak_min(ii))>0, 1)))
                    iMax = 9999;
                else
                    iMax = peak_max(find((peak_max- peak_min(ii))>0, 1));
                end
                new_max = [new_max iMax];
            end
            peak_max = new_max;
            if isempty(peak_min); peak_min = []; end
%         elseif length(peak_max)>length(peak_min)
%             peak_max(end) = [];
%         end
        peak_sel = 0 <= (peak_max-peak_min) & (peak_max-peak_min) <= peak_bw;
        peak_min(~peak_sel) = []; peak_max(~peak_sel) = [];
        
        % find midpoint
        peak_width = 2*[peak_max-peak_min];
        peak_end = peak_min + peak_width; peak_end(peak_end>length(fit_poly)) = length(fit_poly);
        f_mid = (fit_poly(peak_min)+fit_poly(peak_end))*0.5;
        
        [sort_peaks, sort_inds] = sort( fit_poly(peak_max)-f_mid );
        sort_inds(sort_peaks<0) = []; sort_peaks(sort_peaks<0) = [];
        max_peaks = sort_peaks(end-sign(length(sort_peaks)-1):end);
        max_ind = sort_inds(end-sign(length(sort_peaks)-1):end);
        fit_line_x = f_in([peak_max(max_ind); peak_max(max_ind)]);
        fit_line_y = [f_mid(max_ind); fit_poly(peak_max(max_ind))];
        fit_poly = polyval(p, log10(f_in));
        if isempty(max_peaks); max_peaks = 0; end;
        if isempty(fit_line_x); fit_line_x = 0; end;
        if isempty(fit_line_y); fit_line_y = 0; end;
        % output
        fit_params = [];
        fit_params.base_bias = base_bias;
        fit_params.base_exp = base_exp;
        fit_params.peaks = max_peaks;
        fit_params.freqs = fit_line_x(1, :);
        fit_params.order = bestOrder;
        fit_params.p = bestP;
        
        fit_line.fit_line = fit_poly;
        fit_line.xval = fit_line_x;
        fit_line.yval = fit_line_y;
        
%                 figure, plot( log10(x_in)), hold on
%                 plot( fit_poly, 'r');hold on
%                 plot(f_in([peak_min peak_max peak_end]), fit_line([peak_min peak_max peak_end]), 'k+');
        %         plot(fit_line_x, fit_line_y)
end