

function [f, output] = Draw_2A(dt, stimAmp, stimTime,...
                            Imax,Diagram_y,...
                            var, var_names,...
                            cmap1, cmap2, fn)

f = figure;

%   top pane - draw stimulus
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
title(horzcat(fn,': ','Desenzitization in DRG neurons using poke - repoke (Hao 2010, Fig.2A)'));
grid on;


output.stimulus.t=x;
output.stimulus.amp=y;
output.stimulus.defAmp = stimAmp(:,2);

%   middle pane - draw model responses
Ix = nan(1,1);
    
s(2) = axes('OuterPosition', [0 0.4 1 0.4]); 

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

%   bottom pane - draw intensity - response curves
s(3) = axes ('OuterPosition', [0 0 1 0.3]); 

%   boltzmann fit of the Hao2010 2B data
x2a = [0.1:0.1:9];
fig2Boff = Diagram_y - min(Diagram_y);
fig2Bnorm = fig2Boff./max(fig2Boff);
ff = @Boltzmann;
% [param, ~, ~, ~]=fminsearch(ff,[3.6,-2],[],output.stimulus.defAmp,output.model.peakRePoke./min(output.model.peakRePoke));
[paramExp, ~, ~, ~]=fminsearch(ff,[4.8,-1.2],[],output.stimulus.defAmp,fig2Bnorm);
fit2B = 1./(1+exp((paramExp(1)-x2a)/paramExp(2)))+min(Diagram_y)./max(Diagram_y);
fit2BB = fit2B./max(fit2B);
y50Exp = 1./(1+exp((paramExp(1)-paramExp(1))/paramExp(2)))+min(Diagram_y)./max(Diagram_y);
y50Exp = y50Exp./max(fit2B);

%   boltzmann fit of the model 2B/2A model data
fig2modoff = output.model.peakRePoke./min(output.model.peakRePoke);
fig2modoff = fig2modoff - min(fig2modoff);
fig2modnorm = fig2modoff./max(fig2modoff);
[paramMod, ~, ~, ~]=fminsearch(ff,[3,-4],[],output.stimulus.defAmp,fig2modnorm);
% fit2mod = (min(-output.model.peakRePoke)./max(-output.model.peakRePoke))./(1+exp((paramMod(1)-x2a)/paramMod(2)))+min(-output.model.peakRePoke);
fit2mod = range(fig2modoff)./(1+exp((paramMod(1)-x2a)/paramMod(2)))+min(output.model.peakRePoke./min(output.model.peakRePoke));

hold on;
%   plot Boltzmann fit
p5_2B = plot (x2a,fit2BB, '-','LineWidth',2, 'Color',cmap1(7,:));
p6_2B = plot (x2a,fit2mod, '-','LineWidth',2, 'Color',cmap1(2,:));
%   plot x50 values and their projections to the x-axis
plot(paramExp(1),y50Exp,'k^');
line ([paramExp(1), paramExp(1)],[0 y50Exp],'LineStyle','--','Color',cmap1(14,:));
plot(paramMod(1),(min(-output.model.peakRePoke)./max(-output.model.peakRePoke))*0.5...
    +min(-output.model.peakRePoke),'b^');
line ([paramMod(1), paramMod(1)],...
    [0 ((min(-output.model.peakRePoke)./max(-output.model.peakRePoke))*0.5...
    +min(-output.model.peakRePoke))],'LineStyle','--','Color',cmap1(1,:));
hold off;

output.experiment.x50k50 = paramExp;
output.model.x50k50 = paramMod;

hold on;
p1_2B = plot (output.stimulus.defAmp, output.model.peakPoke./min(output.model.peakPoke), '.','Color',cmap1(1,:),'MarkerSize',20);
p2_2B = plot (output.stimulus.defAmp, output.model.peakRePoke./min(output.model.peakRePoke),'*r', 'LineWidth',2,'MarkerSize',10);
p3_2B = plot (output.stimulus.defAmp,Imax./min(Imax), '.','MarkerEdgeColor',cmap1(14,:),'LineWidth',2, 'MarkerSize',20);
p4_2B = plot (output.stimulus.defAmp,Diagram_y./max(Diagram_y),'s','MarkerEdgeColor',cmap1(14,:),'LineWidth',2,'MarkerSize',5); 

text(6,0.1,['exp x_{50} = ', num2str(paramExp(1)),...
    ' exp k = ', num2str(paramExp(2))],...
        'HorizontalAlignment','left','VerticalAlignment','top','Color',cmap1(14,:));

text(6,0.2,['mod x_{50} = ', num2str(paramMod(1)),...
    ' mod k = ', num2str(paramMod(2))],...
        'HorizontalAlignment','left','VerticalAlignment','top','Color',cmap1(2,:));



set(gca,'XTick',output.stimulus.defAmp);
hold off;
ylabel('I/I_{max}');
xlabel('Conditioning stimulus [\mum]');
grid on;


ylim([0 1.01]);
legend ([p1_2B, p3_2B, p2_2B,p4_2B,p5_2B,p6_2B],'cond. peak current model', 'cond. peak current experiment',...
    'model channel availability','experiment channel availability','Boltzmann fit of the experiment','Boltzmann fit of the model','Location','west');



