%%%%%
%   Janez Pre?ern, 2014
%   
%   Function requires:
%       model   ..        model type
%       tT      ..        stimulus time
%       ampT    ..        stimulus amplitude at times tT
%       ProtocolIndex ..  which of the traces/points shall we fit
%       cost6.....        cost functions
%       variables  ..     variable values as inserted by fminsearch
%       variables_names ..names of variables
%       peaksMeasured_recovery .. experimentally measured recovery peaks
%       cw1...cw2..       cost weights
%       wFig6E ..         weights for individual point in the
%                           instensity/response curve

%   Function outputs:
%       outputs.Fig6.model.peakRecovery..  maximum current g at various stimuli amplitudes
%       outputs.c1,.c2,          ...       computed costs

function [out,varargout] = computeFig6 (model,tT,ampT,ProtocolIndex,...
                               cost6Rec_amp,cost6Cost_amp, ...
                                variables,variables_names,dt,...
                                peaksMeasured_recovery,...
                                 cw1,cw2,...
                                 wFig6E)

peakRecovery = nan(1, length(ProtocolIndex));                             
r3 = zeros(1, length(ProtocolIndex));                            
c1 = 0;
c2 = 0;



%%% Fitting the figure 6 - trapezoid stimuli with recovery
for ww = 1:length(ProtocolIndex)            
% parfor ww = 1:length(ProtocolIndex)
    
    %%% computes the g elicited by stimulus
    g = Modeling_DRG_TCM_Engine_v3(model,tT(ww,:),ampT(ww,:),...
        variables,variables_names,dt);
    %%% extract the peaks
    peakRecovery(ww) = min(g(864:end))/min(g(1:864)); % selecting the part of the modeled 
    
    %%%%%%%%%%% C12 compare current traces %%%%%%%%%%%%%%%%%
    if cw1 ~= 0
        e2 = ((g - cost6Rec_amp(ww,:)).^2).*cost6Cost_amp (ww,:);
        e2(isinf(e2)) = NaN;
        r3(ww) = nansum(e2)/(length(e2)-sum(isnan(e2)));
    end;

end;
c1 = nanmean(r3);

%%%%%%%%%%% C3 compare recovery curves %%%%%%%%%%%%%%%%%
%%%         compares the measured peaks against the model
if cw2 ~=0
    c2 = sum((wFig6E(ProtocolIndex)/length(wFig6E(ProtocolIndex))).*...
    (peakRecovery - peaksMeasured_recovery(ProtocolIndex)).^2)./...
    abs((mean(peaksMeasured_recovery(ProtocolIndex)))*length(peakRecovery));
end;

%%%%%  output section
cost = [c1, c2];
out.peakRecovery = peakRecovery;
varargout = {cost};
end