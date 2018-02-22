function F = fit_func3_loglog(x, P, f, f_1)
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

vr = -(sin(x(2))+1)*0.5*x(6);
F=  P - ((x(1)-nansum([x(6)*f; vr*f_1-vr*log10(40)], 1))...
         + x(3)*x(5)*sqrt(2*pi)*normpdf(f,x(4),x(5))...
         + (x(3)-x(7))*x(5)*sqrt(2*pi)*normpdf(f,x(4)+log10(2)+x(8),x(5)));
%(sin(x(7))+1)*0.5*x(3)