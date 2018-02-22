function F = fit_func3_loglog(x, P, f, f_1, base_exp)

% function for fitting a broadband + gaussian

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

% F= P - (x(1)-p_exp*f + x(2)*.04*sqrt(2*pi)*normpdf(f,x(3),.04));
% F= P - (x(1)-p_exp*f + x(2)*.04*sqrt(2*pi)*normpdf(f,x(3),x(4)));
% F= P - (x(1)-x(5)*f + x(2)*.04*sqrt(2*pi)*normpdf(f,x(3),x(4)));
%if x(6)>x(3), x(6) = x(3); end
%if x(4)+log10(2) > f(end)
%    F= P - (x(1)-x(2)*f + x(3)*x(5)*sqrt(2*pi)*normpdf(f,x(4),x(5)));
%else
%F= P - (x(1)-x(2)*f + x(3)*x(5)*sqrt(2*pi)*normpdf(f,x(4),x(5))...
%                    + sqrt((x(3).^2-x(6)))*x(5)*sqrt(2*pi)*normpdf(f,x(4)+log10(2),x(7)));
% F= P - (x(1)-x(2)*f ...
%      + exp(-f*x(9))...
%      + x(3)*x(5)*sqrt(2*pi)*normpdf(f,x(4),x(5))...
%      + x(6)*x(7)*sqrt(2*pi)*normpdf(f,x(8),x(7)));
% F= P - (x(1)-nansum(x(2)*f+x(9)*f2) + x(3)*x(5)*sqrt(2*pi)*normpdf(f,x(4),x(5))...
%     + x(6)*x(7)*sqrt(2*pi)*normpdf(f,x(8),x(7)));
%keyboard
F=  P - ((x(1)-nansum([x(9)*f; x(2)*f_1-x(2)*log10(40)], 1))...
         + x(3)*x(5)*sqrt(2*pi)*normpdf(f,x(4),x(5))...
         + x(6)*x(8)*sqrt(2*pi)*normpdf(f,x(7),x(8)));
% F=  P - ( fit_linear...
%          + x(1)*x(3)*sqrt(2*pi)*normpdf(f,x(2),x(3))...
%          + x(4)*x(6)*sqrt(2*pi)*normpdf(f,x(5),x(6)));

%end

%end
% note: .04*sqrt(2*pi) gives an amplitude of 1 to the Gaussian for x(2)=1;
