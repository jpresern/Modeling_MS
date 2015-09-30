%%%%%
%   Janez Pre?ern, 2014
%   
%   Function requires:
%       model   ..        model type
%       tT      ..        stimulus time
%       ampT    ..        stimulus amplitude at times tT
%       ProtocolIndex ..  which of the traces/points shall we fit
%       cost14.....        cost functions
%       variables  ..     variable values as inserted by fminsearch
%       variables_names ..names of variables
%       peaksMeasured_recovery .. experimentally measured recovery peaks
%       cw1...cw2..       cost weights
%       wFig2   ..        point weights
%   Function outputs:
%       outputs.Fig6.model.peakRecovery..  maximum current g at various stimuli amplitudes
%       outputs.c1,              ...       computed costs

function [out,varargout] = computeFig2 (model,tT,ampT,ProtocolIndex,...
                                variables,variables_names,dt,...
                                peaksMeasured_rePoke,...
                                wFig2)


%   Preparing variables for parfor loop (slicing et al.)
peaksPredicted = nan(size(ProtocolIndex));
peaksPredicted_Repoke = nan(size(ProtocolIndex));
            
%%% Fitting the figure 2A - trapezoid stimuli in poke - repoke fashion

for w14 = 1:length(ProtocolIndex14)    
% parfor w14 = 1:length(ProtocolIndex)

    %%% computes the g elicited by stimulus
    g = Modeling_DRG_TCM_Engine_v3(model,tT(w14,:),ampT(w14,:),...
        variables,variables_names,dt);
    %%% extracts the peaks
    peaksPredicted(w14) = min(g(1:825));  
    peaksPredicted_rePoke(w14) = min(g(826:end)); % selecting the part of the modeled 
                                                  % response we are interested in                
end;

%%%%%%%%%%% C14 compare recovery curves %%%%%%%%%%%%%%%%%
%%%         compares the measured peaks against the model

peaksRepoke = peaksPredicted_rePoke./min(peaksPredicted_rePoke);
c1 = sum((wFig2(ProtocolIndex)/length(wFig2(ProtocolIndex))).*...% preparing the weights for the estimator C2
    (peaksRepoke - peaksMeasured_rePoke(ProtocolIndex)).^2)./...                           % subtracting modeled and recorded peaks
    abs((mean(peaksMeasured_rePoke(ProtocolIndex)))*length(peaksPredicted));                                 % ?normalizing?


%%%%%  output section
cost = c1;
out.peaksRepoke = peaksPredicted_rePoke;
varargout = {cost};
end