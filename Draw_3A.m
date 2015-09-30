function [fig, output] = Draw_3A(dt, stimAmp, stimTime, expAmp, expTime, var,...
                            varInitial, varLimits, var_names,...
                            weights, recWeights, cmap1, cmap2, fn)


% figNum = 1;
fig = figure;
s(1) = axes('OuterPosition', [0 0.8 1 0.2]);  
set(gca,'XTickLabel',[]);
% set(0,'defaultAxesFontName', 'Arial')
% set(0,'defaultTextFontName', 'Arial')

[x,y] = Modeling_GenerateStimulus([0:dt:max(sum(stimTime,2))],...
    stimTime,stimAmp);

hold on;%     ExpData.Fig3A_Recording.ampsmooth{1,a} = [komponenta1', ExpData.Fig3A_Recording.ampi{a}(262:720)', komponenta2'];
for i = 1: size (x,1);
    plot(x(i,:)',y(i,:)','Color',cmap2(i,:),'LineWidth',2);
end;
hold off;

ylabel('x [\mum]');
title(horzcat(fn,': ','Response of DRG neurons to ramp stimuli (Hao 2010, Fig.3A)'));
grid on;

output.stimulus.t=x;
output.stimulus.amp=y;
output.stimulus.ampMax=max(y,[],2);


%%%%%%%%%%%%%%%%%%%%%%%%%% Calculate and draw the responses %%%%%%%%%%%%%%%%%%%%%%
s(2) = axes('OuterPosition', [0 0.4 1 0.4]); 

hold on;
results = [];
for a = 1 : size(stimTime,1)
    expAmp{a} = expAmp{a}';
    res = [];
    [~,~,res]=Modeling_DRG_TCM_Engine_v3(horzcat('DRG_TCM_Model_mh','_Report'),...
        stimTime(a,:),stimAmp(a,:),...
        var,var_names,dt);
   
    p2 = plot(expTime{a},expAmp{a},'LineWidth',1.2,'Color',cmap1(7,:));%cmap(a,:),);
    p1 = plot(res.t,res.g,'-','LineWidth',2,'Color',cmap2(a,:));
    
    output.model.t(a)={res.t};
    output.model.I(a)={res.g};
    [output.model.Imax(a),ixFig3] = min(res.g(1,250:end)); 
    [pk(a), ixFig3Exp(a)] = min(expAmp{a}(1:end));
    plot ((ixFig3+250).*dt,output.model.Imax(a), '.','Color',cmap1(1,:),'MarkerSize',20); 
    plot ((ixFig3Exp(a)).*dt,pk(a),'.','Color','r','MarkerSize',20);
    
    output.experiment.t(a)={expTime(a)};
    output.experiment.I(a)={expAmp(a)};
    output.experiment.Imax(a) = min(expAmp{a}(250:end));
    results{a} = res;
% pause;
end;
output.results = results;
hold off;

hold on;

text(50,-0.15,var_names,'HorizontalAlignment','left','VerticalAlignment','top');
text(60,-0.15,num2str(varInitial'),...
        'HorizontalAlignment','left','VerticalAlignment','top');
text(75,-0.15,num2str(var'),'HorizontalAlignment','left','Color','r','VerticalAlignment','top');
text(90,-0.15,num2str(varLimits),...
    'HorizontalAlignment','left','VerticalAlignment','top',...
    'Color','g');
text(5,0.2,horzcat('Weights: ',num2str(weights)));
text(5,0.1,horzcat('Rec wgs: ',num2str(recWeights)));
hold off;

ylim ([-1.1 0.1]);
xlabel('t [ms]');
ylabel('I/I_{max}');
grid on;

legend([p1,p2],'model','experiment','Location','southeast');

%%%% fit and measure tau rise and tau decay
opt = optimset('MaxFunEval',10000);

[tauMod,~,~,~] = fminsearch(@expDec,3,opt,expTime{13}(1,1:end-(30+ixFig3Exp(a))),-expAmp{13}(1,1+30+ixFig3Exp(a):end));
fitExp = (1-exp(-expTime{13}(1,1:end-(30+ixFig3Exp(a)))./tauMod))-1;
hold on;
p7 = plot(expTime{13}(1,1+30+ixFig3Exp(a):end), fitExp,'Color',cmap1(14,:),'LineWidth',1.5);
hold off;
text(25,-0.9,horzcat('\tau_{fall} = ',num2str(tauMod)),'Color',cmap1(14,:));

[tauExp,~,~,~] = fminsearch(@expDec,3,opt,res.t(1,1:end-(280+ixFig3)),-res.g(1,1+280+ixFig3:end));
fitMod = (1-exp(-res.t(1,1:end-(280+ixFig3))./tauExp))-1;
hold on;
p8 = plot(res.t(1,1+280+ixFig3:end), fitMod,'r-','LineWidth',1.5);
hold off;
text(25,-0.7,horzcat('\tau_{fall} = ',num2str(tauExp)),'Color','r');

grid on;
xlabel('time [ms]');
ylabel('Stim_{50} (\mum)');

output.model.tau = tauMod;
output.experiment.tau = tauExp;

%%%%%%%%%%%%%%%%%%%% plot the intensity response curves 
xx = linspace(0,9,100);

s(3) = axes ('OuterPosition', [0 0 1 0.3]);
hold on;

set(gca,'XTick',stimAmp(:,2));

%%% plot Bolcman first, so the data points are plotted *over* the lines.
ff = @Bolcman;
[paramExp, ~, ~, ~]=fminsearch(ff,[5,0.5],[],output.stimulus.ampMax,...
    output.experiment.Imax'/min(output.experiment.Imax));
fit = 1./(1+exp((paramExp(1)-xx)/paramExp(2)));
p5 = plot (xx,fit,'LineWidth',2,'Color',cmap1(7,:));
plot(paramExp(1),0.5,'*','MarkerSize',10,'MarkerEdgeColor',cmap1(14,:));
line([paramExp(1) paramExp(1)],[0 0.5],'LineStyle','--','Color',cmap1(14,:));

[paramMod, ~, ~, ~]=fminsearch(ff,[5,0.5],[],output.stimulus.ampMax,...
    output.model.Imax'/min(output.model.Imax));
fit = 1./(1+exp((paramMod(1)-xx)/paramMod(2)));
p6 = plot (xx,fit,'LineWidth',2,'Color',cmap1(2,:));
plot(paramMod(1),0.5,'*','MarkerSize',10,'MarkerEdgeColor',cmap1(1,:));
line([paramMod(1) paramMod(1)],[0 0.5],'LineStyle','--','LineWidth',1,'Color',cmap1(1,:));

output.experiment.x50k50 = paramExp;
output.model.x50k50 = paramMod;

p3 = plot (output.stimulus.ampMax,...
    output.experiment.Imax'/min(output.experiment.Imax),...
    's','MarkerEdgeColor',cmap1(14,:),'LineWidth',2,'MarkerSize',5);
p4 = plot (output.stimulus.ampMax,...
    output.model.Imax'/min(output.model.Imax),...
    '.','MarkerEdgeColor',cmap1(1,:),'LineWidth',1,'MarkerSize',20); 

text(6,0.4,['exp x_{50} = ', num2str(paramExp(1)),...
    ' exp k = ', num2str(paramExp(2))],...
        'HorizontalAlignment','left','VerticalAlignment','top','Color',cmap1(14,:));

text(6,0.5,['mod x_{50} = ', num2str(paramMod(1)),...
    ' mod k = ', num2str(paramMod(2))],...
        'HorizontalAlignment','left','VerticalAlignment','top','Color',cmap1(2,:));


hold off;
xlabel ('Stimulus amplitude [\mum]');
ylabel ('I/I_{max}');
grid on;
legend ([p3, p4, p5, p6], 'experiment peak current','model peak current',...
    'Boltzmann fit of experiment','Boltzmann fit of the model','Location','northwest');


% ss{length(f)} = s; s = [];
