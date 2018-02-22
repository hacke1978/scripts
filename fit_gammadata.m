function [fit_params, fit_line] = fit_gammadata(f, f_use4fit, data_base, data_fit)
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
f_in_1=f_in; f_in_1(f_in_1<40) = nan;

if ~isequal(length(f_in), length(x_base)); x_base = x_base'; end
if ~isequal(length(f_in), length(x_in)); x_in = x_in'; end
% fit exponent to base
p=polyfit(log10(f_in),log10(x_base),1);
base_exp  = -p(1);
base_bias =  p(2);

my_options=optimset('Display','off','Algorithm','trust-region-reflective');

x = [];
rms = inf;
finit = linspace(30,80,5);
sdinit = linspace(log10(1./sqrt(2/3)), log10(1./sqrt(19/20)),3);
%slopeinit = linspace(0, -0.5*base_exp,3);
slopeinit = [0 0.25*pi 0.5*pi];
for iInit = 1:length(finit)
    for jSd_init = 1%:length(sdinit)
        for kSlope_init = 1:length(slopeinit)
            [x_tmp, ~, residual]=lsqnonlin(@(x) fit_func3_loglog(x, log10(x_in), log10(f_in), log10(f_in_1)),...
                [ base_bias  slopeinit(kSlope_init) 0   log10(finit(iInit)) sdinit(jSd_init)      base_exp 0 log10(1) ],...
                [ -Inf        -pi                   0   log10(30)           log10(1./sqrt(19/20))  0    0  log10(0.95)  ],...
                [ +Inf        pi                    Inf log10(80)           log10(1./sqrt(2/3))     Inf  Inf log10(1.1)  ],...
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
%x(8) = x(4)+x(8);
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
fit_line = fit_linear + fit_gauss_1 + fit_gauss_2;
