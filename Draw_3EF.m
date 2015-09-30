function [f, output] = Draw_3EF(stimAmp,peakRecovery,...
                                x50k50,...
                                cmap1, cmap2, marmap)

       
% Drawing Fig 3EF: Boltzmann fit

f = figure;

%   conditioning amplitudes
conAmp = [0.0,0.8,1.5,2.2,2.9,3.6,4.3,5.0,5.7,6.4,7.1,7.8,8.5];
%   time points (delays)
delayy = [4.0, 20.0, 40.0, 100.0, 400.0, 1000.0];
Imax = peakRecovery (:,1,1);
Imax = Imax';

%   x values for drawing the fit
xx = linspace(0,9,100);
xxx = linspace(1,1000,1000);

%   conditioning stimuli
x50 = nan (size(stimAmp,2),6);
k50 = nan (size(stimAmp,2),6);

%   offsets to zero and normalizes the input values
% offset = peakRecovery-repmat (max(peakRecovery),[size(peakRecovery,1),1,1]);
offset = peakRecovery;
normal = offset./repmat(min(offset),[size(peakRecovery,1),1,1]);
% normal = -peakRecovery;

%   counts different conditioning amplitudes first
for aa = 1:size(stimAmp,2)

    subplot (3,3,aa);
    x = stimAmp(:,aa);
    yy = normal(:,:,aa);

    x = x(aa:end);
    yy = yy (aa:end,:);
    %   calls the fit and plots for each delay in the matrix
    opt = optimset('MaxFunEval',10000);

    ff = @Bolcman;
    hold on;
    for bb = 1 : size(yy,2);
        y = (yy(:,bb));
        %%% plot points
        
%%%%%   computing x50 & k50 using fminsearch

        [param, ~, ~, ~]=fminsearch(ff,[5,0.5],opt,x(2:end),y(2:end));
              %   puts Stim50 and k @ stim50 into separate variables
        x50(aa,bb) = param(1);
        k50(aa,bb) = param(2);
        fit = 1./(1+exp((param(1)-xx)/param(2)));
        %%% fit
        plot(xx,fit, '-','Color',cmap2(bb,:),'LineWidth',2);
        
        xlabel ('Stimulus [\mum]');
        ylabel ('I/I_{max}');
        xlim ([0 9]);
        ylim ([0 1]);
        title (horzcat ('Cond. stim. amp. ',mat2str(conAmp(aa))));
        grid on;
        
        %%% plot fits
 
        [param, ~, ~, ~]=fminsearch(ff,[5,0.5],opt,stimAmp(:,1),...
                        Imax'/min(Imax));
        fit = 1./(1+exp((param(1)-xx)/param(2)));
        plot(xx,fit, '-','Color',cmap1(7,:),'LineWidth',2);
        plot (x50(aa,bb),0.5, 'Color',cmap1(5,:),'MarkerSize',10,'Marker','x','LineWidth',2);
        plot (x50k50(1),0.5,'Color',cmap1(5,:),'MarkerSize',10,'Marker','x','LineWidth',2);
%         legend('Original','Fit','Location','northwest');
        %%% plot points
        plot(x,y,marmap(bb),'Color',cmap2(bb,:),'LineWidth',2);
        plot (stimAmp(:,1),...
            Imax'/min(Imax),...
            '.','MarkerSize',20,'Color',cmap1(14,:));
    end; 
    hold off;
    set(gca,'XLim',[0 max(stimAmp(:,1))]);
    set(gca,'XTick', conAmp);
    set(gca,'XTickLabel', conAmp);
end;

%%%  Drawing Hao Fig 3EFGHI: Draws the x50 (um) amplitudes against conditioning time (ms)

%   calculate mean values for each conditioning time
%   calculation does not include first line (no slope) and last three lines

%   inserting x50 from control (fig 3A)
x50 = horzcat(repmat(x50k50(1),size(stimAmp,2),1),x50);
% % x50 = x50;
x50(x50<0) = nan;
% %   making sure we do not take stray NaNs in
% x50mean = nanmean(x50(2:9,:),1);

%   offset and normalization
x50min = min(x50,[],2);
x50offset = x50-repmat(x50min,1,7);
x50offsetMax = max(x50offset,[],2);
x50norm = x50offset./repmat(x50offsetMax,1,7);
% x50normMean = nanmean(x50norm(2:9,:),1);

%%%%% PLOT FITS OF ADAPT & INACT
%   open a new figure and plot
f = [f, figure];


s(1) = axes('OuterPosition', [0 0.5 1 0.5]);
set(s(1),'ColorOrder', cmap2(1:size(x50norm,1),:));
fittRange = nan(2, size(x50norm,1));
fitt = nan(size(x50norm,1),size(xxx,2));
tauAct = nan(1,size(x50norm,1));
for dd = 1:size(x50norm,1)
    [tau,~,~,~] = fminsearch(@expApproach,8,opt,[0.00000 delayy],x50norm(dd,:));
    tauAct(dd)=tau;
    fitt(dd,:) = x50offsetMax(dd,1)*(1-exp(-xxx./tau))+x50min(dd,1);
    fittRange(:,dd) = [min(fitt(dd,:)),max(fitt(dd,:))];   
end
semilogx ((repmat([0 delayy],9,1))',x50','xr',...
            (repmat(xxx,9,1))',fitt','-k','LineWidth',2);   
grid on;
xlabel('time [ms]');
ylabel('Stim_{50} (\mum)');
title ('Hao 3F: time course of adaptation shift')

adaptShift = diff(fittRange',[],2);
tauAct(1) = nan;

%%% Drawing 3EFGHI: preparing data Hao Fig 3H - INACTIVATION

%   exctracts peak current values from 
satCur = peakRecovery(size(peakRecovery,1),:,:);
satCur = reshape (satCur,[6, size(peakRecovery,3)])';
%   inserting control current for 9th stimulus amplitude (Fig 3A)
satCur = horzcat(repmat(Imax(13),size (satCur,1),1),satCur);
maxSatCur = max(satCur,[],2);
%   offset to zero
offSatCur = satCur - repmat(maxSatCur,1,7);
%   normalization
normSatCur = offSatCur./repmat(min(offSatCur,[],2),1,7);

s(2) = axes('OuterPosition', [0.0 0.0 1 0.5]);
tauInact = nan (1,size(normSatCur,1));
fittt = nan(size(normSatCur,1),size(xxx,2));
for ff = 2:size(normSatCur,1)
    [tau2,~,~,~] = fminsearch(@expDec,4,opt,[0 delayy],normSatCur(ff,:));
    tauInact(ff) = tau2;
    fittt(ff,:) = exp(-xxx./tau2);    
end;

semilogx ((repmat([0 delayy],9,1))',normSatCur','ob',...
            (repmat(xxx,9,1))',fittt','-k','LineWidth',2); 

grid on;
xlabel('time [ms]');
ylabel('channel availability');
title ('Hao 3F: time course of inactivation')


output.model.tauInact = tauInact;
output.model.tauAct = tauAct;
output.model.adaptShift = adaptShift;
output.model.fittRange = fittRange;

