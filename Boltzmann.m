%   Boltzmann equation is written as in Hao2010
%   Janez Prešern, maj 2014

%%  Fitting the Boltzmann

%   use this if you are using fminsearch
function ret=Bolcman(a,x,y)
    Imax = 1./(1+exp((a(1)-x)/(a(2))));
    ret=sum((y(:)-Imax(:)).^2);

% %   use this if you are using fsqcurvefit
% function ret=Bolcman(a,x,y)
%     ret = 1./(1+exp((a(1)-x)/(a(2))));
    