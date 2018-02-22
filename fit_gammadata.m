function [base_bias, base_exp, event_bias, event_exp_2, gauss_amp, gauss_freq, gauss_std, fit_f2] = ...
     fit_gammadata(f, f_use4fit, data_base, data_fit)
% input:
% f: frequencies
% f_use4fit=[35:57 65:115 126:175 186:200];
% data_base: used to fit exp: (1/f^exp) - enter power (not log)
% data_fit: used to fit weights and gaussian - enter power (not log)
%
% output (exp weight_pwr weight_gauss gamma_freq fit_f2)

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

f_sel=ismember(f,f_use4fit);
x_base=data_base(f_sel);
x_in=data_fit(f_sel);
f_in=f(f_sel);
f_in_1=f_in; f_in_1(f_in_1<40) = nan;
% fit exponent to base
p=polyfit(log10(f_in),log10(x_base)',1);
base_exp  = -p(1);
base_bias =  p(2);

% maxPeak = max(log10(x_in)-event_bias+event_exp*log10(f_in'));
my_options=optimset('Display','off','Algorithm','trust-region-reflective'); % trust-region-reflective

x = [];
rms = inf;
finit = linspace(30,80,5);
sdinit = linspace(log10(1./sqrt(2/3)), log10(1./sqrt(19/20)),3);
slopeinit = linspace(0, -0.5*base_exp,3);
for iInit = 1:length(finit)
    for jSd_init = 1%:length(sdinit)
        for kSlope_init = 1:length(slopeinit)
            [x_tmp, ~, residual]=lsqnonlin(@(x) fit_func3_loglog(x, log10(x_in), log10(f_in'), log10(f_in_1'), base_exp),...
                 [ base_bias  slopeinit(kSlope_init)        0   log10(finit(iInit))   sdinit(jSd_init)       0     log10(finit(iInit)*2) log10(1./sqrt(3/4))    base_exp],... %log10(100)
                 [ 0          -Inf                          0    log10(30)             log10(1./sqrt(19/20))  0     log10(60)             log10(1./sqrt(19/20))  0],... %log10(70)
                 [ Inf        0                             Inf   log10(80)             log10(1./sqrt(2/3))    Inf   log10(f_in(end))      log10(1./sqrt(2/3))   Inf],... %log10(160)
                  my_options);
            if sum(residual.^2)<rms
                rms = sum(residual.^2);
                x = x_tmp;
            end
        end
    end
end
event_bias=x(1);
event_exp_2=x(2);
gauss_amp=x(3);
gauss_freq=x(4);
gauss_std=x(5);
gauss_amp_2=x(6);
gauss_freq_2=x(7);
gauss_std_2=x(8);

% get fits
f_1 = f; f_1(f_1<40) = nan;
fit_linear = event_bias-nansum([x(9)*log10(f) event_exp_2*log10(f_1)-event_exp_2*log10(40)], 2);
fit_gauss_1 = gauss_amp*gauss_std*sqrt(2*pi)*normpdf(log10(f),gauss_freq, gauss_std);
fit_gauss_2 = gauss_amp_2*gauss_std_2*sqrt(2*pi)*normpdf(log10(f),gauss_freq_2, gauss_std_2);
fit_f2 = fit_linear + fit_gauss_1 + fit_gauss_2;
