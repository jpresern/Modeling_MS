%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Janez Presern, Ales Skorjanc, Tomaz Rodic, Jan Benda 2011-2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%   Fitting the exponential decay 
%   Utility function
%   Input
%       tau ..  suggested tau
%       t ..    vector with x values
%       y ..    vector with y values

%   Ouput
%       yy ..   squared difference
function ret=expDec(tau,t,y)
    yy = exp(-t./tau(1));
    ret=sum((y(:)-yy(:)).^2);