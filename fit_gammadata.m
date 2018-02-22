<<<<<<< HEAD
function [base_exp, base_bias, gauss_amp, gauss_freq, fit_f2,gauss_fit] = ...
                    fit_gammadata(f, f_use4fit, data_base, data_fit)

% function fits broadband + gaussian
% [out_exp,w_pwr,w_gauss,gauss_f,fit_f2] = ...
%     fit_gammadata(f,f_use4fit,data_base,data_fit);
%
=======
function [base_bias, base_exp, event_bias, event_exp_2, gauss_amp, gauss_freq, gauss_std, fit_f2] = ...
     fit_gammadata(f, f_use4fit, data_base, data_fit)
>>>>>>> e0d30dd8aa5fcce3792f825ed5367b0eb80756d1
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
<<<<<<< HEAD

if ~isequal(size(f_in),size(x_base))
    x_base=x_base';
end
    
% fit exponent to base
p=polyfit(log10(f_in),log10(x_base),1);

base_exp  = -p(1);
base_bias = p(2);
=======
f_in_1=f_in; f_in_1(f_in_1<40) = nan;
% fit exponent to base
p=polyfit(log10(f_in),log10(x_base)',1);
base_exp  = -p(1);
base_bias =  p(2);
>>>>>>> e0d30dd8aa5fcce3792f825ed5367b0eb80756d1

if ~isequal(size(x_in),size(f_in))
    f_in=f_in';
end

% maxPeak = max(log10(x_in)-event_bias+event_exp*log10(f_in'));
my_options=optimset('Display','off','Algorithm','trust-region-reflective'); % trust-region-reflective
<<<<<<< HEAD
[x]=lsqnonlin(@(x) fit_func3_loglog(x, log10(x_in), log10(f_in)),...
    [ base_bias  base_exp    0   log10(50)   log10(1./sqrt(3/4))    0      log10(1./sqrt(3/4))    ],... %log10(100) 
    [-Inf         -Inf       0   log10(35)   log10(1./sqrt(19/20))  0      log10(1./sqrt(19/20))  ],... %log10(70)    
    [ Inf        base_exp    Inf log10(80)   log10(1./sqrt(2/3))    Inf    log10(1./sqrt(2/3))    ],... %log10(160)   
    my_options);

event_bias  = x(1);
event_exp   = x(2);
gauss_amp   = x(3);
gauss_freq  = x(4);
gauss_std   = x(5);
gauss_amp_2 = x(6);
gauss_std_2 = x(7);

% maxPeak = 1.2 * max(abs(log10(x_in) - (event_bias - event_slope*log10(f_in'))));
% [x]=lsqnonlin(@(x) fit_func3_loglog(x,log10(x_in),log10(f_in'),base_exp),...
%     [ 0     0       log10(50) log10(1.5)   base_exp  0        log10(gauss_freq)   log10(1.5)  ],... %
%     [-Inf   0       log10(35) log10(0)    -Inf       0        log10(gauss_freq)   log10(0)    ],... %
%     [ Inf   maxPeak log10(80) log10(3.5)   Inf       maxPeak  log10(gauss_freq)   log10(2.5)  ],... %
%     my_options);

% second fit
% my_options=optimset('Display','off','Algorithm','trust-region-reflective');
% [x]=lsqnonlin(@(x) fit_func3_loglog(x,log10(x_in),log10(f_in'),out_exp),...
%     [0 0 log10(50)],[-Inf 0 log10(35)],[Inf Inf log10(80)],...
%     my_options);

% fit to data in log-space
% fit_f2=w_pwr-out_exp*log10(f) + ...
%     w_gauss*.04*sqrt(2*pi)*normpdf((f),gauss_f,.04);
% fit_f2=w_pwr-out_exp*log10(f) + ...
%     w_gauss*.04*sqrt(2*pi)*normpdf(log10(f),gauss_f, gauss_w);
% fit_f2=w_pwr-gauss_s*log10(f) + ...
%     w_gauss*.04*sqrt(2*pi)*normpdf(log10(f),gauss_f, gauss_w);
% figure
fit_linear = event_bias-event_exp*log10(f);
fit_gauss_1 = gauss_amp*.04*sqrt(2*pi)*normpdf(log10(f),gauss_freq, gauss_std);
fit_gauss_2 = gauss_amp_2*.04*sqrt(2*pi)*normpdf(log10(f),gauss_freq+log10(2), gauss_std_2);
fit_f2 = fit_linear + fit_gauss_1 + fit_gauss_2; 

%AP I added this for more convenient saving/plotting
gauss_fit.amp   = gauss_amp;
gauss_fit.freq  = 10.^gauss_freq;
gauss_fit.std   = gauss_std;
gauss_fit.amp2  = gauss_amp_2;
gauss_fit.std2  = gauss_std_2;
gauss_fit.freq2 = 2*(10.^gauss_freq);
gauss_fit.event_exp   = event_exp;
gauss_fit.event_bias  = event_bias;
gauss_fit.base_exp    = base_exp;
gauss_fit.base_bias   = base_bias;
gauss_fit.fit_f2      = fit_f2;
% loglog(f, 10.^(fit_f2), 'color', [1 0 0 ]/2)
% hold on
% loglog(f, 10.^(event_bias-event_slope*log10(f)), 'color', [1 0 0 ]/2)
% hold on
% fit_f2=...
%     ...% event_bias-event_slope*log10(f) + ...
%     gauss_amp*.04*sqrt(2*pi)*normpdf(log10(f),gauss_freq, gauss_std) + ...
%     gauss_amp_2*.04*sqrt(2*pi)*normpdf(log10(f),gauss_freq+log10(2), gauss_std_2);
% loglog(f, 10.^(fit_f2), 'color', [1 0 0 ]/2)
=======

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
>>>>>>> e0d30dd8aa5fcce3792f825ed5367b0eb80756d1
