%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Janez Presern, Ales Skorjanc, Tomaz Rodic, Jan Benda 2011-2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [f] = Draw_3CD (dt, stim,...
                            ampMax, Imax, x50k50,...
                            var, var_names,...
                            cmap1, cmap2, marmap, fn)



%%%% Conditioning stimuli loop
for l = 5%1:9        %   conditioning higher than 5.7 was not done as it is not possible to fit with boltzmann (not enough points)
    

f = figure;

set (gca, 'Visible', 'off');

%%%%%%%%%%%% drawing stimulus
s(1) = axes ('OuterPosition', [0 0.8 1 0.2]);
set(gca,'XTickLabel',[]); 

stimAmp = [];

hold on;
for val1 = 1:size (stim(1,l).t,1)
    for val3 = 1:size (stim(1,l).t,3)
        x0 = 0;
        y0 = 0;
        for val2 = 1:7
            x1 = x0 + stim(1,l).t(val1,val2,val3);
            y1 = y0 + stim(1,l).amp(val1,val2,val3);
            line([x0 x1],[y0 y1], 'Color',cmap2(val3,:));
            x0 = x1;
            y0 = y1;
        end;
        stimAmp(val1, val3) = sum (stim(1,l).amp(val1,1:4,val3));
    end;
end;
xlim([0,1100]);
ylabel ('stimulus amplitude [\mum]');
hold off;
grid on;
title(horzcat(fn,': ','Desensitization process (Hao 2010, Fig.3CDE)'));


%%%%%%%%%%%%% drawing model responses
s(2) = axes ('OuterPosition', [0 0.4 1 0.4]);

%   prepares the variables for the both the conditioning stimulus peak and
%   peak due to the test stimulus
peakInitial = NaN(size(stim(1,l).t,1),size(stim(1,l).t,3));
peakRecovery = peakInitial;
ix1 = peakInitial;
% ix2 = peakInitial;
ix2 = [];


hold on;

%   prepare variable to store current, computed in the response to
%   conditioning as a reference for peak extraction.
condResp = [];    

%%%%% Amplitude loop - start at slikca value as it mathces the line of NaN
%%%%% Delay loop    
for n = 1 : size (stim(1,l).t,3) 
    
%%%%% amplitudes
    for m = l : size(stim(1,l).t,1)
    
  
        %   clears the variables, makes sure no left-over artefacts from previous run are drawn (or
        %   extracted)
        results.g = [];
        results.t = [];
        %    checks if the data imput contains a number or NaN
        if ~isnan(stim(1,l).t(m,1:3,n))
            [~,~,results]=Modeling_DRG_TCM_Engine(horzcat('DRG_TCM_Model_mh','_Report'),...
                stim(1,l).t(m,:,n),stim(1,l).amp(m,:,n),...
                var,var_names,dt);
            
        %   calculates index to grab peaks in the correct time window
            
            
            %   stores the current computed in the response to the
            %   conditioning stimuli, storing only max. delay current to ensure the
            %   proper length of the vector for further computation;
            %   However it stores peak values and index

            if m == l
                condResp = results.g;
                [peakInitial(m,n), ix1(m,n)] = min(condResp);
                %%% problem with no peaks at the 0 stimuli
%                     if dummy3(m,n) == length(condResp)
                if ix1(m,n) >= int64(sum(stim(1,l).t(m,1:3,n))/dt);
                    ix1(m,n) = int64(sum(stim(1,l).t(m,1:3,n))/dt);
                end
                
            else
                %   if no value is stored, copy the value along the
                %   dimensions of the matrix for simpler use
                peakInitial(m,:) = peakInitial(l,:);
                ix1(m,:) = ix1(l,:);
            end;
            
            %   recovery peak extraction
            [pks, ix2] = findpeaks (-results.g(ix1(m,n):end),'npeaks',1);
            ix2 = ix1(m,n) + ix2;
        end;
        peakRecovery (m,n) = 0;
        if ~isempty(pks)
%             peakRecovery (m,n) = -pks-condResp(ix2+ix1(m,n));
            peakRecovery (m,n) = -pks-condResp(ix2);
        end;
        plot(results.t,results.g,'-','Color',cmap2(n,:));
        plot(results.t(ix2),-pks,'.b','MarkerSize',15);

        
    end;    
end;
xlim([0,1100]);
xlabel ('time [ms]');
ylabel ('I/I_{max}');
grid on;

hold off;


%%%%%%%%%% drawing intensity response curves Hao 3E

s(3) = axes ('OuterPosition', [0 0 0.5 0.3]);

conAmp = [0.0,0.8,1.5,2.2,2.9,3.6,4.3,5.0,5.7,6.4,7.1,7.8,8.5];
delayy = [4.0, 20.0, 40.0, 100.0, 400.0, 1000.0];
xx = linspace(0,9,100);
xxx = linspace(1,1000,1000);

%   conditioning stimuli
x50 = nan (size(stimAmp,2),6);
k50 = nan (size(stimAmp,2),6);

%   offsets to zero and normalizes the input values
% offset = peakRecovery-repmat (max(peakRecovery),[13,1,1]);
offset = peakRecovery;
normal = offset./repmat(min(offset),[13,1,1]);

markers  = ['o'];%,'+','*','x','^','v'];
hold on;
plot (ampMax,Imax'/min(Imax),'s','Color',cmap1(14,:),'MarkerSize',5, 'LineWidth',2);
for i = 1:size(stimAmp,2) 
    plot ((stimAmp(:,i)), (peakRecovery(:,i)./min(Imax)), marmap(i),'Color',cmap2(i,:),'LineWidth',2);
end;
set(gca,'XTick',conAmp);
hold off;

xlabel ('Stimulus [\mum]');
ylabel ('I/I_{max}');
grid on;
ylim ([0 1.01]);
xlim ([0 9]);

s(4) = axes ('OuterPosition', [0.5 0.0 0.5 0.3]);

for l = 5%1:size(stimAmp,2)
%     fig(1,l) = figure;
    
    x = stimAmp(:,l);
%     yy = normal(:,:,l);
    yy = normal;

    x = x(l:end);
    yy = yy (l:end,:);
    %   calls the fit and plots for each delay in the matrix
    opt = optimset('MaxFunEval',10000);
%     ret = @(a,x)1./(1+exp((a(1)-x)/(a(2))));
    ff = @Boltzmann;
    hold on;
    for bb = 1 : size(yy,2);
        y = (yy(:,bb));
        %%% plot points
        
%%%%%   computing x50 & k50 using fminsearch

        [param, ~, ~, ~]=fminsearch(ff,[5,0.5],opt,x(2:end),y(2:end));
%         %   computing x50 & k50 using lsqcurvefit 
%         [param,~,~,~,~] = lsqcurvefit (ret,[5;0.5],x(2:end),y(2:end));
        
        %   puts Stim50 and k @ stim50 into separate variables
        x50(l,bb) = param(1);
        k50(l,bb) = param(2);
        fit = 1./(1+exp((param(1)-xx)/param(2)));
        %%% fit
        plot(xx,fit, '-','Color',cmap2(bb,:),'LineWidth',2);
        
        xlabel ('Stimulus [\mum]');
        ylabel ('I/I_{max}');
        xlim ([0 9]);
        ylim ([0 1]);
%         title (horzcat ('Cond. stim. amp. ',mat2str(conAmp(l)),' \mum Delay ', mat2str(delayy(bb))));
        grid on;
        
        %%% plot fits
 
        [param, ~, ~, ~]=fminsearch(ff,[5,0.5],opt,ampMax,...
                        Imax'/min(Imax));
        fit = 1./(1+exp((param(1)-xx)/param(2)));
        plot(xx,fit, '-','Color',cmap1(7,:),'LineWidth',2);
        plot (x50(l,bb),0.5, 'Color',cmap1(5,:),'MarkerSize',10,'Marker','x','LineWidth',2);
        plot (x50k50(1),0.5,'Color',cmap1(5,:),'MarkerSize',10,'Marker','x','LineWidth',2);
%         legend('Original','Fit','Location','northwest');
        %%% plot points
        plot(x,y,marmap(bb),'Color',cmap2(bb,:),'LineWidth',2);
        plot (ampMax,...
            Imax'/min(Imax),...
            's','MarkerSize',5,'Color',cmap1(14,:),'LineWidth',2);
    end; 
    hold off;
    set(gca,'XTick', conAmp);
    set(gca,'XTickLabel', conAmp);
end;




ss{length(f)} = s; s = [];

end;