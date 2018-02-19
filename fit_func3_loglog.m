function F = fit_func3_loglog(x,P,f)

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
F= P - (x(1)-x(2)*f + x(3)*.04*sqrt(2*pi)*normpdf(f,x(4),x(5))...
                    + x(6)*.04*sqrt(2*pi)*normpdf(f,x(4)+log10(2),x(7)));
% note: .04*sqrt(2*pi) gives an amplitude of 1 to the Gaussian for x(2)=1;
