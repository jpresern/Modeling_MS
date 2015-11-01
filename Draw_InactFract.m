%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Janez Presern, Ales Skorjanc, Tomaz Rodic, Jan Benda 2011-2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%   Function draws the inactivated fraction in responses to desensitization 
%   with a prepulse (like Hao & Delmas, Fig 2A). 

%   Function requires:
%       dt..                sampling rate
%       stimamp..           stimulus description (amplitudes changes)
%       stimTime..          stimulus description (time durations of amplitude changes

%       var          ..     variable values as inserted by fminsearch
%       var_names           ..names of variables
%       cmap1, cmap2 ..     color maps
%       fn ..               filename (Results_XX) to identify the model
%                           iteration
%       
%   Function outputs:
%       f             ..    figure handle
%       output.stimulus.t ..generated time course of stimulus
%       output.stimulus.amp.generated amplitude course of stimulus
%       output.stimulus.defAmp..colected maximum amplitudes of stimulus
%       output.model.t ..   time line of model response
%       output.model.I ..   time course of current(amplitude) model
%                           response
%       output.experiment.x50k50 .. mid point and slope of experimentally 
%                           obtained I-R curve for inactivation (Hao Fig2)
%       output.model.x50k50.mid point and slope for model generated I-R
%                           curve for inactivation
%       output.model.peakPoke .. model peak current of the prepulse
%       output.model.peakRePoke..model peak current of the test stimulus

function [f, output] = Draw_InactFract(dt, stimAmp, stimTime,...
                            var, var_names,...
                            cmap1, cmap2, fn)

f = figure;

%%%%%%%%%%%%%%%%%%%%%%%%%% Calculate and draw the stimulus %%%%%%%%%%%%%%%%%%%%%%
s(1) = axes('OuterPosition', [0 0.8 1 0.2]);
set(gca,'XTickLabel',[]);

[x,y] = Modeling_GenerateStimulus([0:dt:max(sum(stimTime,2))],...
    stimTime,stimAmp);

hold on;
for i = 1: size (x,1)
    plot(x(i,:)',y(i,:)','Color',cmap2(i,:),'LineWidth',2);  
end;
hold off;

ylabel('x [\mum]');
xlim ([0 90]);
title(horzcat(fn,': ','Inactivated fraction with double-pulse paradigm'),'interpreter','none');
grid on;

output.stimulus.t=x;
output.stimulus.amp=y;
output.stimulus.defAmp = stimAmp(:,2);

%%%%%%%%%%%%%%%%%%%%%%%%%% Calculate and draw the responses %%%%%%%%%%%%%%%%%%%%%%    
s(2) = axes('OuterPosition', [0 0.35 1 0.4]); 

hold on;
for a = 1 : size(stimTime,1)
    results.t = [];
    results.g = [];
    [~,~,results]=Modeling_DRG_TCM_Engine(horzcat('DRG_TCM_Model_mh','_Report'),...
        stimTime(a,:),stimAmp(a,:),...
        var,var_names,dt);
    
    p1 = plot(results.t,results.g,'-','LineWidth',2,'Color',cmap2(a,:));
    
    output.model.t(a)={results.t};
    output.model.I(a)={results.g};
   
    %   harvesting the peak values
    firstPeak = int64(sum(stimTime(a,1:3))/dt);
    secondPeak = int64(sum(stimTime(a,1:4))/dt);

    %   grabs the peaks
    [output.model.peakPoke(a), ixPoke] = min(results.g(1:firstPeak));
    [output.model.peakRePoke(a), IxRe] = min(results.g(secondPeak:end));
    
    plot(ixPoke*dt,output.model.peakPoke(a),'.b','MarkerSize',20);
    plot((length(1:secondPeak)+IxRe).*dt,output.model.peakRePoke(a),'*r','MarkerSize',7);

end;
hold off;
ylim ([-1.1 0]);
xlim ([0 90]);
xlabel('t [ms]');
ylabel('I/I_{max}');
grid on;

%%%%%%%%%%%%%%%%%%%% Calculate and draw the intensity response curves 
s(3) = axes ('OuterPosition', [0 0 1 0.3]); 

step = [0.1:0.1:9];
ff = @Boltzmann;

%   boltzmann fit of the model chan. inactivation and availability data
fig2modoff = output.model.peakRePoke./min(output.model.peakRePoke);
fig2modoff = fig2modoff - min(fig2modoff);
fig2modnorm = fig2modoff./max(fig2modoff);
[paramMod, ~, ~, ~]=fminsearch(ff,[3,-4],[],output.stimulus.defAmp,fig2modnorm);
fit2mod = range(fig2modoff)./(1+exp((paramMod(1)-step)/paramMod(2)))+min(output.model.peakRePoke./min(output.model.peakRePoke));

hold on;
%   plot Boltzmann fit
p6_2B = plot (step,fit2mod, '-','LineWidth',2, 'Color',cmap1(2,:));
%   plot x50 values and their projections to the x-axis
% plot(paramMod(1),(min(-output.model.peakRePoke)./max(-output.model.peakRePoke))*0.5...
%     +min(-output.model.peakRePoke),'b^');
% line ([paramMod(1), paramMod(1)],...
%     [0 ((min(-output.model.peakRePoke)./max(-output.model.peakRePoke))*0.5...
%     +min(-output.model.peakRePoke))],'LineStyle','--','Color',cmap1(1,:));
hold off;

output.model.x50k50 = paramMod;

hold on;
p1_2B = plot (output.stimulus.defAmp, output.model.peakPoke./min(output.model.peakPoke), '.','Color',cmap1(1,:),'MarkerSize',20);
p2_2B = plot (output.stimulus.defAmp, output.model.peakRePoke./min(output.model.peakRePoke),'*r', 'LineWidth',2,'MarkerSize',10);

text(6,0.2,['mod x_{50} = ', num2str(paramMod(1)),...
    ' mod k = ', num2str(paramMod(2))],...
        'HorizontalAlignment','left','VerticalAlignment','top','Color',cmap1(2,:));

set(gca,'XTick',output.stimulus.defAmp);
hold off;
ylabel('I/I_{max}');
xlabel('Conditioning stimulus [\mum]');
grid on;

ylim([0 1.01]);

legend ([p1_2B, p2_2B,p6_2B],'cond. peak current model',...
    'model channel availability','Boltzmann fit of the model','Location','west');


