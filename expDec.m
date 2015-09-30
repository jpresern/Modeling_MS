%   Exponential decay
%   Janez Prešern, maj 2014

%%  Fitting the exponential approach

function ret=expDec(tau,t,y)
    yy = exp(-t./tau(1));
    ret=sum((y(:)-yy(:)).^2);