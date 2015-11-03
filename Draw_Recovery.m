%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Janez Presern, Ales Skorjanc, Tomaz Rodic, Jan Benda 2011-2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%   Function draws the recovery from desensitization

%   Function requires:
%       dt..                sampling rate
%       stimamp..           stimulus description (amplitudes changes)
%       stimTime..          stimulus description (time durations of amplitude changes
%       var  ..             fitted variables and constants
%       var_names ..        names of variables
%       cmap1,cmap2 ..      color maps
%       fn ..               filename (Results_XX) to identify the model
%                           iteration

%   Function outputs:
%       f             ..    figure handle
%   output.model.peakInitial..peak currents to the prepulse
%   output.model.peakRecovery..peak currents to the test stimulus


function [f, output] = Draw_Recovery (dt, stimAmp, stimTime,...
                                    var,var_names,...
                                    cmap1, cmap2, fn)

f = figure; 

%%%%%%%%%%%%%%%%%%%%% Draw the stimulus %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
s(1) = axes('OuterPosition', [0 0.8 1 0.2]);  
set(gca,'XTickLabel',[]);
[x,y] = Modeling_GenerateStimulus([0:dt:max(sum(stimTime,2))],...
    stimTime,stimAmp);

hold on;
for i = 1: size (x,1);
    plot(x(i,:)',y(i,:)','Color',cmap2(i,:),'LineWidth',2);
end;
hold off;

ylabel('x [\mum]');
title(horzcat(fn,': ','Recovery of MS channels from inactivation'),'interpreter','none');
grid on;
xlim([0 120]);

%%%%%%%%%%%%%%%%%%% Draw the model responses %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
s(2) = axes('OuterPosition', [0 0.3 1 0.5]);
peakInitial = nan(1,8);
peakRecovery = peakInitial;

hold on;
for  a = 1 : size(stimTime,1)
    
    [~,~,results]=Modeling_DRG_TCM_Engine(horzcat('DRG_TCM_Model_mh','_Report'),...
        stimTime(a,:),stimAmp(a,:),...
        var,var_names,dt);
    peakInitial(a) = min(results.g(1:1001));
    [peakRecovery(a), ixRec] = min(results.g(1001:end));
    
    p1 = plot(results.t,results.g,'-','LineWidth',2,'Color',cmap2(a,:));
    
    output.model.t(a)={results.t};
    output.model.I(a)={results.g};
    plot((ixRec+1001).*dt,peakRecovery(a), '*r', 'MarkerSize',7);

end;
ylim ([-1.1 0]);
xlim([0 120]);
hold off;

output.model.peakInitial = peakInitial;
output.model.peakRecovery = peakRecovery;

xlabel('t [ms]');
ylabel('I/I_{max}[%]');
grid on;

ss{length(f)} = s; s = [];

%%%%%%% draw intensity/response curve %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

s(3) = axes('OuterPosition', [0 0 1 0.3]);

hold on;

p3 = plot(stimTime(1:end,5)',peakRecovery./min(peakInitial),'*r','LineWidth',2,'MarkerSize',10);

hold off;

axis([0 10 0 1.1]);
xlabel('stimulus delay [ms]');
ylabel('I/I_{max}[%]');
grid on;
xlim ([0 45]);

%   fit and plot the exponential approach on the collected and normalized
%   modeled peaks
timRec = 0.1:0.1:40;

offsetPeakRecovery = (-peakRecovery-min(-peakRecovery));
normPeakRecovery = offsetPeakRecovery./max(offsetPeakRecovery);
opt = optimset('MaxFunEval',10000);
tauRec = [];
[tauRec,~,~,~] = fminsearch(@expApproach,8,opt,stimTime(1:end,5)',normPeakRecovery);
hold on;
fitRec = range(peakRecovery).*(1-exp(-timRec./tauRec))+(1-max(offsetPeakRecovery));
p1 = plot(timRec,fitRec, '-','Color',cmap1(2,:),'LineWidth',2);
hold off;

hold on;
text(15,0.65,['\tau = ',num2str(tauRec)],'HorizontalAlignment','left','Color',cmap1(2,:));
hold off;

legend ([p3, p1],'model','Fit of the model','Location','southeast');



