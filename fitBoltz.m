%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Janez Pre?ern, 2014
%   Function to be applied over the matrix row-/columnwise instead of using
%   for loop. Hopefully it will cut down time for calculations


%   output is normally into a single layer matrix, therefore k50 output
%   is disabled by default

% function [x50,k50] = fitBolc (x,y)
function x50 = fitBoltz(x,y)
    ff = @Boltzmann;
    opt = optimset('MaxFunEval',10000);
    [param, ~] = fminsearch(ff,[5,0.5],opt,x(1:end),y(1:end));
    x50 = param(1);
%     k50 = param(2);
    