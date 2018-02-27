function F = fit_loglog_exp(x, P, f, flagStep, expParam)
% fit a broadband exp + 2 gaussians
%     Cem Uran - cem.uran@esi-frankfurt.de
%
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
%

if flagStep
    if nargin == 4
        F=  P - ((x(1)-x(2)*exp(-f*x(3))));
    elseif nargin == 5
        F=  P - ((expParam(1)-expParam(2)*exp(-f*expParam(3)))...
            + x(1)*x(3)*sqrt(2*pi)*normpdf(f,x(2),x(3))...
            + x(4)*x(3)*sqrt(2*pi)*normpdf(f,x(2)+log10(2)+x(5),x(3)));
    end
else
    F=  P - ((x(1)-x(2)*exp(-f*x(3)))...
        + x(4)*x(6)*sqrt(2*pi)*normpdf(f,x(5),x(6))...
        + x(7)*x(9)*sqrt(2*pi)*normpdf(f,x(8),x(9)));
end
