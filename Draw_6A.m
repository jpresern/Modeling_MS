%%  Drawing the Fig6
stimAmp = [0.0,0.8,1.5,2.2,2.9,3.6,4.3,5.0,5.7,6.4,7.1,7.8,8.5];
dur = [5, 10, 20, 50, 100, 200, 500, 600, 1000];
xx = linspace(0,9,100);
xxx = linspace(1,1100,1100);

% 3x3 DRG responses like Fig6A for all stimulus amplitudes

f = [f, figure];

s(2) = axes('OuterPosition', [0 0.4 1 0.4]); 

aa = 1:1:13;
bb = 1:1:9;
output.Fig6A.model.peakPoke = nan(length(aa),length(bb));
output.Fig6A.model.peakPoke = nan(length(aa),length(bb));
% matlabpool local 3
hold on;
for b = bb
    subplot (3,3,b);
    hold on;
    for a = aa
    
        %   computes
        [~,~,results]=Modeling_DRG_TCM_Engine_v3(horzcat(InitialFitParam.Model,'_Report'),...
            ExpData.Fig6A_Stimulus.t(a,:,b),ExpData.Fig6A_Stimulus.amp(a,:,b),...
            var,var_names,InitialFitParam.dt);
        %   plots
        p1 = plot(results.t,results.g,'-','Color',cmap(a,:),'LineWidth',2);
        %   grabs first (control) peak and second test peak
        firstPeak = int64(sum(ExpData.Fig6A_Stimulus.t(a,1:3,b))/InitialFitParam.dt);
        secondPeak = int64(sum(ExpData.Fig6A_Stimulus.t(a,1:4,b))/InitialFitParam.dt);
        output.Fig6A.model.peakPoke(a,b) = min(results.g(1:firstPeak));
        [output.Fig6A.model.peakRePoke(a,b), Ix] = min(results.g(secondPeak:end));

    %     line ([output.Fig6A.model.tIMax(b) output.Fig6A.model.tIMax(b)],[output.Fig6A.model.Imax(b) 0],'Color','r');
        
    end;
    hold off;
    set(gca,'Ylim',[-1.1 0]);
%     xlim([0 1050]);
    xlabel ('time [ms]');
    ylabel ('I/I_{max}');
%     pause;
end;

% plot (output.Fig6A.model.tIMax, output.Fig6A.model.Imax, '*r','MarkerSize',5);
% matlabpool close;

%  Fig 6B2

%   prepares empty storage space - ten units (control + 9 tests)
x50 = nan (10,1);
k50 = nan (10,1);

%   prepares function
ff = @Bolcman;

%   opens new figure
f = [f, figure];
hold on;

%   calculates & draws control
y = output.Fig6A.model.peakPoke(:,1)./repmat(min(output.Fig6A.model.peakPoke(:,1)),13,1);
opt = optimset('MaxFunEval',10000);
plot(stimAmp,y,'o','Color','k','Marker','o','MarkerSize',5);
[param, ~, ~, ~]=fminsearch(ff,[5,0.5],opt,stimAmp,y);
fit = 1./(1+exp((param(1)-xx)/param(2)));
plot(xx,fit, '-','Color','k','LineWidth',1);

%   puts Stim50 and k @ stim50 into separate variables
x50(1,1) = param(1);
k50(1,1) = param(2);
plot (x50(1,1),0.5,'Color','k','Marker','x','MarkerSize',14);

%   calculates & plots tests
for aa = 1:size(output.Fig6A.model.peakPoke,2)
    x = stimAmp(:);
    y = output.Fig6A.model.peakRePoke(:,aa)./repmat(min(output.Fig6A.model.peakRePoke(:,aa)),13,1);
    %   calls the fit and plots for each delay in the matrix
    opt = optimset('MaxFunEval',10000);
    plot(x,y,'o','Color',cmap(aa,:));
    [param, ~, ~, ~]=fminsearch(ff,[5,0.5],opt,x,y);
    %   puts Stim50 and k @ stim50 into separate variables
    x50(aa+1,1) = param(1);
    k50(aa+1,1) = param(2);
    %   draws calculated fit
    fit = 1./(1+exp((param(1)-xx)/param(2)));
    plot(xx,fit, '-','Color',cmap(aa,:),'LineWidth',1);
    plot (x50(aa+1,1),0.5,'Color',cmap(aa,:),'Marker','x','MarkerSize',14);
end;
hold off;
xlabel ('Stimulus [\mum]');
ylabel ('I/I_{max}');
xlim ([0 9]);
ylim ([0 1]);
grid on;

set(gca,'XTick', stimAmp);
set(gca,'XTickLabel', stimAmp);

%  Fig 6C

x50offset = x50 - min(x50);
x50norm = x50offset./max(x50offset);

f = [f, figure];
hold on;
plot ([0 dur],x50, '-or');
[tau,~,~,~] = fminsearch(@expApproach,8,opt,[0 dur],x50norm);
fitt = max(x50offset)*(1-exp(-xxx./tau))+min(x50);
plot(xxx,fitt, '-','Color','k','LineWidth',1);
hold off;
xlim ([0 1100]);
xlabel ('time of conditioning [ms]');
ylabel ('Stim50 [\mum]');
title('Hao2010: Figure 6C');
tau
