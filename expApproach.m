%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Janez Presern, Ales Skorjanc, Tomaz Rodic, Jan Benda 2011-2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Fitting the exponential approach
%   Utility function
%   Input
%       tau ..  suggested tau
%       t ..    vector with x values
%       y ..    vector with y values

%   Ouput
%       yy ..   squared difference

function ret=expApproach(tau,t,y)
    yy = 1-exp(-t./tau);
    ret=sum((y(:)-yy(:)).^2);