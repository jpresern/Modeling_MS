%   Exponential approach
%   Janez Pre�ern, maj 2014

%%  Fitting the exponential approach

function ret=expApproach(tau,t,y)
    yy = 1-exp(-t./tau);
%     xxx = linspace (0,1000,1000);
%     fitt = 1-exp(-xxx./tau);
%     hold on;
%     plot(xxx,fitt, '-','LineWidth',1);
%     hold off;
    ret=sum((y(:)-yy(:)).^2);