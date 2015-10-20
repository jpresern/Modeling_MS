%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Janez Presern, Ales Skorjanc, Tomaz Rodic, Jan Benda 2011-2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Boltzmann equation used for fitting is written as in Hao2010

%   use this if you are using fminsearch
function ret=Boltzmann(a,x,y)
    Imax = 1./(1+exp((a(1)-x)/(a(2))));
    ret=sum((y(:)-Imax(:)).^2);

% %   use this if you are using fsqcurvefit
% function ret=Bolcman(a,x,y)
%     ret = 1./(1+exp((a(1)-x)/(a(2))));
    