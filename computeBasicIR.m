%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Janez Presern, Ales Skorjanc, Tomaz Rodic, Jan Benda 2011-2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Function computes the main stimulus (poke-and-hold). It compares the
%   model response to the published experiment traces. At the same time it
%   also compares the I-R curve of the model to the publisehd one.

%   Function requires:
%       model   ..        model type
%       tT      ..        stimulus time
%       ampT    ..        stimulus amplitude at times tT
%       ProtocolIndex ..  which of the traces/points shall we fit
%       cost1..           cost functions
%       variables  ..     variable values as inserted by fminsearch
%       variables_names ..names of variables
%       peaksMeasured   ..peaks measured experimentally
%       cw1...cw5..       cost weights
%       wT  ..         weights for individual trace
%       wP ..         weights for individual point in the
%                           instensity/response curve

%   Function outputs:
%       output.basicIR.model.iMax..  maximum current g at various stimuli amplitudes
%       output.basicIR.model.amp ..  amplitude used to get iMax
%       output.c1,.c2,.c3...       computed costs

function [out,varargout] = computeBasicIR (model,tT,ampT,ProtocolIndex,...
                               cost1Rec_amp,cost1Cost_amp, cost1Rec_t,...
                                variables,variables_names,dt,...
                                peaksMeasured,...
                                 cw1,cw2,cw3,cw4,cw5,...
                                 wT,wP)

%   variable prealocation
iMax = nan(1,length(ProtocolIndex));
r = zeros(1, length(ProtocolIndex));
cS = 0;
cSS = 0;
c2 = 0;
c5 = 0;

%   loop computing model response to stimuli, extracting the features and
%   comparing them to the recordings
% for w = 1:length(ProtocolIndex)
parfor w = 1:length(ProtocolIndex)
    
    g = Modeling_DRG_TCM_Engine(model,tT(w,:),ampT(w,:),...
        variables,variables_names,dt);
    %%% extracts the peaks
    iMax(w) = min(g(1,250:end)); % not from teh begining to avoid bumps
    
    %%%%%%%%%%% C1 compare current traces %%%%%%%%%%%%%%%%%                
    % Cost is calculated by subtracting actual recording from modeled and 
    % multiplied by cost function defined % separately
    if cw1 ~= 0
        e2 = ((g - cost1Rec_amp(w,:)).^2).*cost1Cost_amp(w,:); 
        e2(isinf(e2)) = NaN;                                           
        r(w) = nansum(e2)/(length(e2)-sum(isnan(e2))); 
    end
    %%%%%%%%%%%%% C7 restrict current traces from basic I-R to approach 0 mV in steady-state %%%%%%%%%%%%%%%%
    if cw3 ~= 0
        indS = dsearchn(cost1Rec_t(w,:)',70);
        if g(indS) < -0.05
            cS(w) = 1000;
        end;
    end;

    %%%%%%%%%%% C9 keep current traces before the stimulus in basic I-R above 0.02 %%%%%%%%%%%%
    if cw4 ~= 0
        indSS = dsearchn(cost1Rec_t(w,:)',18);
        if g(indSS) > -0.01
            cSS(w) = 1000;
        end;
    end;
end;
            
c1 = nanmean(r.*wT(1:length(r)));   % calculates penalties over all recordings
c3 = nanmean(cS);
c4 = nanmean(cSS);


%%%%%%%%%%%%%%%% C2 compare onset response curves %%%%%%%%%%%%%%%%%%%%%%%%
%%% comparing current peaks from the model as produced by protocol from
%%% basic I-R
if cw2 ~= 0
    c2 = sum((wP(ProtocolIndex)/length(wP(ProtocolIndex))).*...% preparing the weights for the estimator C2
        (iMax - peaksMeasured(ProtocolIndex)).^2)./...               % subtracting modeled and recorded peaks
        abs((mean(peaksMeasured(ProtocolIndex)))*length(iMax));       % ?normalizing?
end;

%%%%%%%%%%%%%%%% C10 penalize if min peak goes below -1.1 or if it doesn't go below -0.9  %%%%%%%%%%%%%%%%
if cw5 ~= 0
    if min(iMax) < -1.1 || iMax(end-1)/iMax(end) < 0.95
        c5 = 1000;
    end;
    if min(iMax) > 0.9
        c5 = 1000;
    end;
end;

%%%%%% output section
if strcmp(model,'DRG_TCM_Model_mh_Report')
    out.t = res.t;
    out.g = res.g;
end;    
cost = [c1, c2, c3, c4, c5];
out.Imax = iMax;
varargout = {cost};
end      
