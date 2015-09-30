%%%%%
%   Janez Presern, 2014
%   
%   Function requires:
%       model   ..        model type
%       stimulus      ..  stimulus - structure containing 3D matrix with 
%                           stim times & stimulus amplitudes
%       variables  ..     variable values as inserted by fminsearch
%       variables_names ..names of variables
%       control_amps ...  amplitude from control stimulus ExpData.Fig3A....
%       control_pks ..    peaks from the ExpData.Fig3A.modeled.iMax
%       cw1...cw2..       cost weights
%       wFig2   ..        point weights
%   Function outputs:
%       outputs.Fig6.model.peakRecovery..  maximum current g at various stimuli amplitudes
%       outputs.c1,              ...       computed costs ...


function [out, varargout] = compute3GHI (model,stimulus,...
                                varr,varr_names,dt,...
                                control_pks,control_amps,...
                                tauFig3GH, Fig3I,...
                                cw1,cw2,cw3,Fig3I_PW)
                            

                            
%   some definition
delayy = [4.0, 20.0, 40.0, 100.0, 400.0, 1000.0];
% xxx = linspace(1,1000,1000);
%   enter the control peak values and computes Boltzmann fit over it
opt = optimset('MaxFunEval',10000);
ff = @Bolcman;
[par, ~, ~, ~]=fminsearch(ff,[5,0.5],opt,control_amps,...
                     control_pks/min(control_pks));
x50_cont = par(1);
k50_cont = par(2);

respReduction_measured = Fig3I(:,3);
adaptShift_measured = Fig3I(:,3);

% wFig3I_inact = wFig3I (:,2);
% wFig3I_adapt = wFig3I (:,3);

tauFig3GH = reshape(tauFig3GH(:,2:3),1,8);
wFig3I_adapt = Fig3I_PW (:,2);
wFig3I_inact = Fig3I_PW (:,1);

%%%%%%%%%%%%% calculating the model responses

%   prepares the variables for the both the conditioning stimulus peak and
%   test stimulus    
peakInitial = NaN(size(stimulus(1,1).t,1),size(stimulus(1,1).t,3),size(stimulus,2));
peakRecovery = peakInitial;
x50 = NaN(size(stimulus,2),size(stimulus(1,1).t,3));
k50 = x50;
%%%%%% Conditioning stimulus loop
for l = 1 : size(stimulus,2)
% parfor l = 1 : size(stimulus,2)
%   It turned out to be very efficient to have the parfor loop as the
%   outernmost for loop. The calculation time was cut down from 300-500 s
%   to ~ 150 s including opening the matlab pool (10 s) and closing it (5
%   s)
    %   parfor: prepare dummy variables for peakInitial and peakRecovery
    dummy1=nan(size(stimulus(1,l).t,1),size(stimulus(1,l).t,3));
    dummy2=nan(size(stimulus(1,l).t,1),size(stimulus(1,l).t,3));
    %   parfor: prepare dummy variables for index ix1 & ix2
    dummy3=nan(size(stimulus(1,l).t,1),size(stimulus(1,l).t,3));
    condResp = [];
%     opt = optimset('MaxFunEval',10000);
%     ret = @(a,x)1./(1+exp((a(1)-x)/(a(2))));
%    l
    
    %%%%%% Delay loop
    for n = 1 : size (stimulus(1,l).t,3) 
%         n
        %%%%%% Amplitude loop, starts at value l to skip the rows with NaN    
        for m = l : size(stimulus(1,l).t,1)  
%             m
            %   clears the variables, makes sure no artefacts are drawn (or
            %   extracted)            
            ix2 = [];
            pks = [];
            %    checks if stimulus is a number or NaN
            if ~isnan(stimulus(1,l).t(m,1:3,n)) 
                %    calculates the model responses
                g = Modeling_DRG_TCM_Engine_v3(model,stimulus(1,l).t(m,:,n),...
                                stimulus(1,l).amp(m,:,n),varr,varr_names,dt);

            %   calculates index to grab peaks in the correct time window
                %   stores the current computed in the response to the
                %   conditioning stimuli, storing only max. delay current to ensure the
                %   proper length of the vector for further computation;
                %   However it stores peak values and index
               
                if m == l
                    condResp = g;
                    [dummy1(m,n), dummy3(m,n)] = min(condResp);

                    %%% problem with no peaks at the 0 stimuli
%                     if dummy3(m,n) == length(condResp)
                    if dummy3(m,n) >= int64(sum(stimulus(1,l).t(m,1:3,n))/dt);
                        dummy3(m,n) = int64(sum(stimulus(1,l).t(m,1:3,n))/dt);
                    end;
                else
                    %   if no value is stored, copy the value along the
                    %   dimensions of the matrix for simpler use
                    dummy1(m,:) = dummy1(l,:);
                    dummy3(m,:) = dummy3(l,:);
                end;
                

                %   recovery peak extraction
                [pks, ix2] = findpeaks (-g(dummy3(m,n):int64(sum(stimulus(1,l).t(m,1:end-2,n))/dt)),...
                                'npeaks',1);
                ix2 = dummy3(m,n) + ix2;
            end;
            dummy2(m,n) = 0;
            %   subtracting the value of the conditioning response (offset)
            if ~isempty(pks)
%                 dummy2(m,n) = -pks-condResp(ix2+dummy3(m,n));
                dummy2(m,n) = -pks-condResp(ix2);
            end;
%         plot(g);
%         hold on;
%         plot(ix2,-pks,'*r');
%         plot(dummy3(m,n),dummy1(m,n),'xb');
%         hold off;
        end;
    end;
    peakInitial(:,:,l) = dummy1;
    peakRecovery(:,:,l) = dummy2;
end;

%%%%%%%%%%% calculates Boltzmann when done with the peakRecovery
%%%%%%%%%%% calculations - it is faster in the separated loop


normal = peakRecovery./repmat(min(peakRecovery),[13,1,1]);
x50 = NaN(size(stimulus,2),size(stimulus(1,1).t,3));
k50 = x50;
for l = 1:size(normal,3)
% parfor l = 1:size(normal,3)
    x = control_amps(l:end)';
    x50_temp = nan(size(normal,2),1);
    k50_temp = x50_temp;
    for n = 1 : size (stimulus(1,l).t,3) 
        [param, ~, ~, ~]=fminsearch(ff,[5,0.5],opt,x(1:end),normal(l:end,n,l));
        x50_temp(n,1) = param(1);
        k50_temp(n,1) = param(2);
    end
    x50(l,:) = x50_temp;
    k50(l,:) = k50_temp;
end



%%%%%%%%%%% preparing Hao 3G: computation of Tau activation and adaptive
%%%%%%%%%%% shift

%   inserting the control 
%   THE CONTROL IS SMALLER THAN THE SMALLEST PEAK... USE THAT ONE!
x50 = horzcat(repmat(x50(1,1),size(stimulus,2),1),x50);
% x50 = horzcat(repmat(x50_cont,size(stimulus,2),1),x50);
%   offset and normalization
x50(x50<0) = nan;
x50min = min(x50,[],2);
x50offset = x50-repmat(x50min,1,7);
x50offsetMax = max(x50offset,[],2);
x50norm = x50offset./repmat(x50offsetMax,1,7);
% x50normMean = nanmean(x50norm(2:9,:),1);

%   fit and plot the exponential approach on the collected and normalized stim50 
% fittRange = nan(2, size(x50norm,1));
opt = optimset('MaxFunEval',10000);%,'Display','iter');
tauAct = nan(1,size(x50norm,1));
% hold on;
for dd = 2:size(x50norm,1)
% parfor dd = 2:size(x50norm,1)
    [tau,~,~,~] = fminsearch(@expApproach,8,opt,[0.00000 delayy],x50norm(dd,:));
    %%%% JUST TAKE THE AVERAGE OF THE LAST THREE POINTS INSTEAAAAAD to compute the channel availabilty/inactivated fraction!
%     fitt = x50offsetMax(dd,1)*(1-exp(-xxx./tau))+x50min(dd,1);
%     plot(xxx,fitt);
    tauAct (dd)=tau;
%     fittRange(:,dd) = [min(fitt),max(fitt)]; 

end;
% hold off;
%   the NaN fit returns attempted value, therefore we correct it.
tauAct(1) = NaN;
%%%% calculating the adaptive shift
% fittRange = diff(fittRange);
adaptShift = mean(x50offset(:,4:end),2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% preparing data Hao Fig 3H and inactviated
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% fraction/channel availability



%   exctracts peak current values from 
satCur = peakRecovery(13,:,:);
satCur = reshape (satCur,[6, size(peakRecovery,3)])';

%   inserting max current from control (Fig 3A)
satCur = horzcat(repmat(control_pks(13),size (satCur,1),1),satCur);

% Is it neccessary to do the offset?
maxSatCur = max(satCur,[],2);
%   offset to zero
offSatCur = satCur - repmat(maxSatCur,1,7);

% %   no offset
% offSatCur = satCur;

%   normalization
normSatCur = offSatCur./repmat(min(offSatCur,[],2),1,7);


%   fit and plot the exponential decay on the collected and normalized peak currents
tauInact = nan (1,size(normSatCur,1));
% hold on;
for ff = 2:size(normSatCur,1)
% parfor ff = 2:size(normSatCur,1)
    [tau2,~,~,~] = fminsearch(@expDec,4,opt,[0 delayy],normSatCur(ff,:));
%     fittt = exp(-xxx./tau2);
%     plot(xxx,fittt);%, '-','Color',cmap(ff,:),'LineWidth',1);
    tauInact (ff)=tau2(1);
end;
% hold off;

%%%%% calculating ration of reduction/availability for 3I


peakSat = peakRecovery(13,6,:);
peakSat = reshape(peakSat,size(peakSat,3),1);
%   calculates the ratio betweem max. current (at 8.5) and max current at
%   8.5 um on the different conditioning stimuli using the max delay
chanAvailability = peakSat./min(peakSat,[],1);
respReduction = 1 - chanAvailability;



%%%%% cost calculation


%   fiting 3GH: program compares distribution of TauInact & TauAdapt from
%   the model with the distribution of the measured Taus. 
%   TauInact & TauAdapt are merged together in a single vector for that
%   very purpose. Should H0 be rejected, there is a steep penalty to pay.


if cw1 ~= 0
    c1 = sum((tauFig3GH/length(tauFig3GH)).*...
        ([tauAct(6:9),tauInact(6:9)] - tauFig3GH)).^2./...
        abs(mean([tauAct(6:9),tauInact(6:9)])*length(tauFig3GH));
    
%     [h,~] = ttest ([tauAct(6:9),tauInact(6:9)],tauFig3GH);
%     if h ~= 0
%         c1 = 1000;
%     else
%         c1 = h;
%     end;
else
    c1 = 0;
end;

%   fiting 3I adaptation
if cw2 ~= 0
    c2 = sum((wFig3I_adapt/length(wFig3I_adapt)).*...
        (adaptShift - adaptShift_measured)).^2./...
        abs(mean(adaptShift_measured)*length(adaptShift));
else
    c2 = 0;
end;

%   fiting 3I inactivation
if cw3 ~= 0
    c3 = sum((wFig3I_inact/length(wFig3I_inact)).*...
        (respReduction - respReduction_measured)).^2./...
        abs(mean(respReduction_measured)*length(respReduction));
else
    c3 = 0;
end;




%%%%%% output section

cost = [c1, c2, c3];
output.Fig3A.iMax = control_pks;
out.peakInitial = peakInitial;
out.peakRecovery = peakRecovery;
out.x50 = x50;
out.k50 = k50;
out.x50_cont = x50_cont;
out.k50_cont = k50_cont;
out.tauAct = tauAct;
out.tauInact = tauInact;
out.fittRange = adaptShift;
out.chanAvailability = chanAvailability;
out.respReduction = respReduction;
varargout = {cost};                            
                            
end