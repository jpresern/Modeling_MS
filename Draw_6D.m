function [f, output] = Draw_6D (dt, stimAmp, stimTime,...
                                    Diagram_x, Diagram_y,...
                                    var,var_names,...
                                    cmap1, cmap2, fn)

f = figure; 

%%%%%%%%%%%%%%%%%%%%% Draw the stimulus
s(1) = axes('OuterPosition', [0 0.8 1 0.2]);  
set(gca,'XTickLabel',[]);
[x,y] = Modeling_GenerateStimulus([0:dt:max(sum(stimTime,2))],...
    stimTime,stimAmp);

hold on;
for i = 1: size (x,1);
    plot(x(i,:)',y(i,:)','Color',cmap2(i,:),'LineWidth',2);
%     plot(x(i,:)',y(i,:)');
end;
hold off;

ylabel('x [\mum]');
title(horzcat(fn,': ','Recovery of MS channels from inactivation (Hao 2010, Fig.6D)'));
grid on;
xlim([0 120]);

%%%%%%%%%%%%%%%%%%% Draw the model responses %%%%%%%%%%%%%%%%
s(2) = axes('OuterPosition', [0 0.3 1 0.5]);
peakInitial = nan(1,8);%nan(1,size(ExpData.Fig6D_10ms_Stimulus.t2,1));
peakRecovery = peakInitial;

hold on;
for  a = 1 : size(stimTime,1)
    
    [~,~,results]=Modeling_DRG_TCM_Engine(horzcat('DRG_TCM_Model_mh','_Report'),...
        stimTime(a,:),stimAmp(a,:),...
        var,var_names,dt);
    peakInitial(a) = min(results.g(1:1001));
    [peakRecovery(a), ixRec] = min(results.g(1001:end));
    
    p1 = plot(results.t,results.g,'-','LineWidth',2,'Color',cmap2(a,:));
%     p2 =  plot(ExpData.Fig6D_10ms_Recording.ti{a},ExpData.Fig6D_10ms_Recording.ampi{a},'--','Color',cmap(a,:),'LineWidth',LineWidth);
    
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
% legend([p1,p2],'model','experiment','Location','southeast');

ss{length(f)} = s; s = [];

%%%%%%% draw intensity/response curve Hao 2010, Fig 6E %%%%%%%%%

s(3) = axes('OuterPosition', [0 0 1 0.3]);

hold on;

p3 = plot(stimTime(1:end,5)',peakRecovery./min(peakInitial),'*r','LineWidth',2,'MarkerSize',10);
p4 = plot(Diagram_x,Diagram_y./max(Diagram_y),...
    's', 'LineWidth',2, 'MarkerSize',5,'Color',cmap1(14,:));

hold off;

axis([0 10 0 1.1]);
xlabel('stimulus delay [ms]');
ylabel('I/I_{max}[%]');
grid on;
xlim ([0 45]);

%   fit and plot the exponential approach on the collected and normalized
%   modeled peaks
timRec = 0.1:0.1:40;
% peakRecovery = -peakRecovery./min(peakRecovery);
offsetPeakRecovery = (-peakRecovery-min(-peakRecovery));
normPeakRecovery = offsetPeakRecovery./max(offsetPeakRecovery);
opt = optimset('MaxFunEval',10000);
tauRec = [];
[tauRec,~,~,~] = fminsearch(@expApproach,8,opt,stimTime(1:end,5)',normPeakRecovery);
hold on;
fitRec = range(peakRecovery).*(1-exp(-timRec./tauRec))+(1-max(offsetPeakRecovery));
p1 = plot(timRec,fitRec, '-','Color',cmap1(2,:),'LineWidth',2);
hold off;

%   fit and plot the exponential approach on the collected and normalized
%   experimental peaks
timExp = 0.1:0.1:40;
offsetPeakExp = (Diagram_y-min(Diagram_y));
normPeakExp = offsetPeakExp./max(offsetPeakExp);
opt = optimset('MaxFunEval',10000);
tauExp = [];
[tauExp,~,~,~] = fminsearch(@expApproach,8,opt,Diagram_x,normPeakExp);
hold on;
fitExp = range(Diagram_y).*(1-exp(-timExp./tauExp))+(1-max(offsetPeakExp));
p2 = plot(timExp,fitExp, '-','Color',cmap1(7,:),'LineWidth',2);
hold off;

hold on;
text(15,0.65,['\tau = ',num2str(tauRec)],'HorizontalAlignment','left','Color',cmap1(2,:));
text(15,0.45,['\tau = ',num2str(tauExp)],'HorizontalAlignment','left','Color',cmap1(14,:));
hold off;

legend ([p4, p2, p3, p1],'experiment','Fit of the experiment','model','Fit of the model','Location','southeast');



