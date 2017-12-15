
%% https://en.wikipedia.org/wiki/Visual_angle
%% V = 2 * arctan(S/2*D)

%% distance to monitor  =  64 cm
%% screensize = [1682 1052] px = [477x298] mm
%% 3.526 is ratio of px to mm
%% radian to degree conversion is done by *57.2958


function [ S ] = degDistances( D)


% D = 2*atan((S/3.526)/)*57.2958;  %% output is in degrees 


S = 3.556*(830*2)*tan(D/(2*57.2958)); %% output is in pixels 
end

