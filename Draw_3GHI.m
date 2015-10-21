%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Janez Presern, Ales Skorjanc, Tomaz Rodic, Jan Benda 2011-2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%   Function computes and draws time constants of inactivation and adaptation
%   Function also draws contributions of adaptation and inactivation to the
%   desensitization.

%   Function requires:
%       exp3GHIInactAdapt.. published time constants of inactivation and
%                           adaptation (Hao 2010, Fig 3)

%       exp3IInactAdapt ..  published contributions of adaptation and
%                           inactivation (Hao 2010, Fig 3I)
%       tauAct  ..          time constant of adaptation obtained from model
%       tauInact ..         time constant of inactivation from model
%       adaptShift ..       adaptive shift from the model
%       peakRecovery ..     peak currents from desentization tests
%       cmap1       ..      color maps

%   Function outputs:
%       f             ..    figure handle
%   output.model.paramInactExp.. mid point and slope for Boltzmann fit of
%                                published experimental contribution of
%                                inactivation (Hao 2010, Fig3I)
%   output.model.paramAdaptExp ..mid point and slope for Boltzmann fit of
%                                published experimental contribution of
%                                adaptation (Hao 2010, Fig3I)
%   output.model.paramInact ..   mid point and slope for Boltzmann fit of
%                                contribution of inactivation from the model
%   output.model.paramAdapt ..   mid point and slope for Boltzmann fit of
%                                contribution of adaptation from the model
%   output.model.ratioInact ..   ratio between inactivated and available
%                                fraction from the model

function [f, output] = Draw_3GHI (exp3GHInactAdapt,exp3IInactAdapt,...
                                    tauAct, tauInact, adaptShift,...
                                    peakRecovery,...
                                    cmap1)


%   conditioning amplitudes
conAmp = [0.0,0.8,1.5,2.2,2.9,3.6,4.3,5.0,5.7,6.4,7.1,7.8,8.5];
Imax = peakRecovery (:,1,1);
Imax = Imax';


f = figure;

s(1) = axes('OuterPosition', [0 0.7 1 0.3]);
hold on;
p3 = plot (exp3GHInactAdapt(:,1),exp3GHInactAdapt(:,2), '.', 'LineWidth',2, 'MarkerSize',20,'Color',cmap1(14,:));
p4 = plot(exp3GHInactAdapt(:,1),exp3GHInactAdapt(:,3),'x','LineWidth',2, 'MarkerSize',10,'Color',cmap1(14,:));
p1 = plot (conAmp(1:size(tauAct,2)), tauAct, 'rx', 'LineWidth',2,'MarkerSize',10);
p2 = plot(conAmp(1:size(tauAct,2)), tauInact, 'b.', 'LineWidth',2,'MarkerSize',20);
hold off;
%text(0.8,0.75,num2str(tauAct));
legend ([p1,p2,p3,p4],'model \tau adaptive shift','model \tau inactivation','Hao2010: \tau adaptive shift','Hao2010: \tau inactivation','Location','southwest');
ylim ([0 10]);
xlim ([3 6]);
ylabel ('\tau [ms]');
xlabel ('Conditioning stimulus [\mum]');
title ('Model adaptive shift & inactivation');
set(gca,'XTick', conAmp);
set(gca,'XTickLabel', conAmp);
grid on;

%%%  Drawing figure 3EFGHI: diagram 3I
%   Inactivation: data needed are from the 3D fit (as it shows the inactivation
%   Adaptive shift: data needed are min & max from 3F for all delays
%                   (values are obtained from the fit)

%   peakInitial is harvested out of Fig3A instead of Fig3CD. Computations
%   for Fig3CD do not neccessary include all the amplitudes.
Imax = Imax';
peakRecovery = peakRecovery(13,6,:);
peakRecovery = reshape(peakRecovery,size(peakRecovery,3),1);
%   calculates the ratio betweem max. current (at 8.5) and max current at
%   8.5 um on the different conditioning stimuli using the max delay
% ratio = peakRecovery./min(Imax,[],1)  
ratio = -(min(Imax,[],1)-peakRecovery)  


s(2) = axes('OuterPosition', [0 0.35 1 0.3]);
[AX, H3, H4] = plotyy (exp3IInactAdapt(:,1), exp3IInactAdapt(:,2), exp3IInactAdapt(:,1), exp3IInactAdapt(:,3));
% [AX, H1, H2] = plotyy (conAmp(1:size(ratio,1)),1-ratio, conAmp(1:size(ratio,1)), adaptShift./max(adaptShift));%adaptShift));
ylim(AX(1),[0, 1]);ylim(AX(2),[0, 1]);% max(adaptShift)]);
xlim(AX(1),[0 6.0]);xlim(AX(2),[0 6.0]);
set(H3, 'Color',cmap1(14,:),'Marker','.','LineStyle','none','MarkerSize',20,'LineWidth',2);
set(H4, 'Color',cmap1(14,:), 'Marker','x','LineStyle','none','MarkerSize',10,'LineWidth',2);
set(AX(1),'Ycolor','blue');
set(AX(2),'Ycolor','red');
xlabel(AX(1),'Conditioning stimuli [\mum]');
ylabel(AX(1),'Inactivated fraction');
ylabel(AX(2),'Normalized adaptive shift');
title(AX(1),'Hao Fig3I: inactivated fraction and adaptive shift as a function of conditioning');
grid (AX(1), 'on');
grid (AX(2), 'on');
    
set(AX(1),'XTick', conAmp);
set(AX(1),'XTickLabel', conAmp);
set(AX(2),'XTick', conAmp);
set(AX(2),'XTickLabel', conAmp);
set(AX(1),'YTick',[0,0.25,0.5,0.75,1.0]);
set(AX(2),'YTick',[0,0.25,0.5,0.75,1.0,1.25,1.5,1.75,2.0,2.5,3.0]);

%   boltzmann fit of the Hao2010 2B data
x2a = [0.1:0.1:9];
hold on;
ff = @Boltzmann;
[paramInactExp, ~, ~, ~]=fminsearch(ff,[4.8,1.2],[],exp3IInactAdapt(:,1),exp3IInactAdapt(:,2));
fit3IinactExp = 1./(1+exp((paramInactExp(1)-x2a)/paramInactExp(2)));
plot(x2a, fit3IinactExp,'Color',cmap1(7,:),'LineWidth',2);

% adaptShiftNorm = adaptShift./max(adaptShift);

[paramAdaptExp, ~, ~, ~]=fminsearch(ff,[4.8,1.2],[],exp3IInactAdapt(:,1),exp3IInactAdapt(:,3));
fit3IadaptExp = 1./(1+exp((paramAdaptExp(1)-x2a)/paramAdaptExp(2)));
plot(x2a, fit3IadaptExp,'Color',cmap1(7,:),'LineWidth',2)


plot(paramInactExp(1),0.5,'^','Color',cmap1(14,:),'LineWidth',2);
line ([paramInactExp(1), paramInactExp(1)],[0 0.5],'LineStyle','--','Color',cmap1(14,:));
plot(paramAdaptExp(1),0.5,'^','Color',cmap1(14,:),'LineWidth',2);
line ([paramAdaptExp(1), paramAdaptExp(1)],[0 0.5],'LineStyle','--','Color',cmap1(14,:));

hold off;

legend([H3, H4], 'Hao2010 inactivated fraction', 'Hao2010 adaptive shift','location','northwest');

%%% Draw activation shift and inactivation
s(3) = axes('OuterPosition', [0 0 1 0.3]); 
[AX, H3, H4] = plotyy (exp3IInactAdapt(:,1), exp3IInactAdapt(:,2), exp3IInactAdapt(:,1), exp3IInactAdapt(:,3));
% [AX, H1, H2] = plotyy (conAmp(1:size(ratio,1)),1-ratio, conAmp(1:size(ratio,1)), adaptShift./max(adaptShift));%adaptShift));
ylim(AX(1),[0, 1]);ylim(AX(2),[0, 1]);% max(adaptShift)]);
xlim(AX(1),[0 6.0]);xlim(AX(2),[0 6.0]);
set(H3, 'Color',cmap1(14,:),'Marker','.','LineStyle','none','MarkerSize',20,'LineWidth',2);
set(H4, 'Color',cmap1(14,:), 'Marker','x','LineStyle','none','MarkerSize',10,'LineWidth',2);
set(AX(1),'Ycolor','blue');
set(AX(2),'Ycolor','red');
xlabel(AX(1),'Conditioning stimuli [\mum]');
ylabel(AX(1),'Inactivated fraction');
ylabel(AX(2),'Normalized adaptive shift');
title(AX(1),'Hao Fig3I: inactivated fraction and adaptive shift as a function of conditioning');
grid (AX(1), 'on');
grid (AX(2), 'on');
    
set(AX(1),'XTick', conAmp);
set(AX(1),'XTickLabel', conAmp);
set(AX(2),'XTick', conAmp);
set(AX(2),'XTickLabel', conAmp);
set(AX(1),'YTick',[0,0.25,0.5,0.75,1.0]);
set(AX(2),'YTick',[0,0.25,0.5,0.75,1.0,1.25,1.5,1.75,2.0,2.5,3.0]);

%   boltzmann fit of the Hao2010 2B data
x2a = [0.1:0.1:9];
hold on;
ff = @Boltzmann;
[paramInactExp, ~, ~, ~]=fminsearch(ff,[4.8,1.2],[],exp3IInactAdapt(:,1),exp3IInactAdapt(:,2))
fit3IinactExp = 1./(1+exp((paramInactExp(1)-x2a)/paramInactExp(2)));
plot(x2a, fit3IinactExp,'Color',cmap1(7,:),'LineWidth',2);

% adaptShiftNorm = adaptShift./max(adaptShift);

[paramAdaptExp, ~, ~, ~]=fminsearch(ff,[4.8,1.2],[],exp3IInactAdapt(:,1),exp3IInactAdapt(:,3))
fit3IadaptExp = 1./(1+exp((paramAdaptExp(1)-x2a)/paramAdaptExp(2)));
plot(x2a, fit3IadaptExp,'Color',cmap1(7,:),'LineWidth',2)


plot(paramInactExp(1),0.5,'^','Color',cmap1(14,:),'LineWidth',2);
line ([paramInactExp(1), paramInactExp(1)],[0 0.5],'LineStyle','--','Color',cmap1(14,:));
plot(paramAdaptExp(1),0.5,'^','Color',cmap1(14,:),'LineWidth',2);
line ([paramAdaptExp(1), paramAdaptExp(1)],[0 0.5],'LineStyle','--','Color',cmap1(14,:));

hold off;

hold on;
%s(2) = axes('OuterPosition', [0 0.5 1 0.5]);
[AX, H1, H2] = plotyy (conAmp(1:size(ratio,1)),ratio, conAmp(1:size(ratio,1)), adaptShift./max(adaptShift));%adaptShift));
ylim(AX(1),[0, 1]);ylim(AX(2),[0, 1]);% max(adaptShift)]);
xlim(AX(1),[0 6.0]);xlim(AX(2),[0 6.0]);
ylim(AX(1),[0, 1]);ylim(AX(2),[0, 1]);% max(adaptShift)]);
xlim(AX(1),[0 6.0]);xlim(AX(2),[0 6.0]);
set(H1, 'Color','blue','Marker','.','LineStyle','none','MarkerSize',20','LineWidth',2);
set(H2, 'Color','red','Marker','x','LineStyle','none','MarkerSize',10,'LineWidth',2);
set(AX(1),'Ycolor','blue');
set(AX(2),'Ycolor','red');
xlabel(AX(1),'Conditioning stimuli [\mum]');
ylabel(AX(1),'Inactivated fraction');
ylabel(AX(2),'Normalized adaptive shift');
title(AX(1),'Hao Fig3I: inactivated fraction and adaptive shift as a function of conditioning');
grid (AX(1), 'on');
grid (AX(2), 'on');
    
set(AX(1),'XTick', conAmp);
set(AX(1),'XTickLabel', conAmp);
set(AX(2),'XTick', conAmp);
set(AX(2),'XTickLabel', conAmp);
set(AX(1),'YTick',[0,0.25,0.5,0.75,1.0]);
set(AX(2),'YTick',[0,0.25,0.5,0.75,1.0,1.25,1.5,1.75,2.0,2.5,3.0]);

%   turning off the axis on the opposite side of tick labels
set(AX(1),'Box','Off');
set(AX(2),'Box','Off');
hold off;
%   Boltzmann fit for the curves

%   boltzmann fit of the Hao2010 2B data
x2a = [0.1:0.1:9];
hold on;
ff = @Boltzmann;
[paramInact, ~, ~, ~]=fminsearch(ff,[4.8,1.2],[],conAmp(1:size(ratio,1)),ratio)
fit3Iinact = 1./(1+exp((paramInact(1)-x2a)/paramInact(2)));
plot(x2a, fit3Iinact,'b','LineWidth',2);

% adaptShiftNorm = adaptShift./max(adaptShift);

[paramAdapt, ~, ~, ~]=fminsearch(ff,[4.8,1.2],[],conAmp(1:size(adaptShift,1)),adaptShift./max(adaptShift))
fit3Iadapt = 1./(1+exp((paramAdapt(1)-x2a)/paramAdapt(2)));
plot(x2a, fit3Iadapt,'r','LineWidth',2);

plot(paramInact(1),0.5,'b^','LineWidth',2);
plot(paramAdapt(1),0.5,'r^','LineWidth',2);
line ([paramInact(1), paramInact(1)],[0 0.5],'LineStyle','--','Color','b');
line ([paramAdapt(1), paramAdapt(1)],[0 0.5],'LineStyle','--','Color','r');
hold off;

legend([H1, H2, H3, H4], 'model inactivated fraction','model adaptive shift','Hao2010 inactivated fraction', 'Hao2010 adaptive shift','location','northwest');

output.model.paramInactExp = paramInactExp;
output.model.paramAdaptExp = paramAdaptExp;
output.model.paramInact = paramInact;
output.model.paramAdapt = paramAdapt;
output.model.ratioInact = ratio;
